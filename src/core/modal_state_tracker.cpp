#include "modal_state_tracker.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace BellOptimizer {
namespace {
constexpr std::size_t kShapeSamples = 32;
constexpr std::size_t kShapeComponents = 3;

bool makeShape(const ModalMode& mode, std::vector<double>& out,
               std::vector<int>& harmonics) {
    const auto n = mode.outer_surface_radial.size();
    if (n < 2 || mode.outer_surface_axial.size() != n) return false;

    if (!mode.outer_surface_harmonics.empty() &&
        n > std::numeric_limits<std::size_t>::max() / mode.outer_surface_harmonics.size())
        return false;
    const std::size_t fourier_size = mode.outer_surface_harmonics.empty()
        ? 0u : n * mode.outer_surface_harmonics.size();
    const bool has_fourier = !mode.outer_surface_harmonics.empty() &&
        mode.outer_surface_fourier_radial.size() == fourier_size &&
        mode.outer_surface_fourier_azimuthal.size() == fourier_size &&
        mode.outer_surface_fourier_axial.size() == fourier_size;
    harmonics = has_fourier ? mode.outer_surface_harmonics : std::vector<int>{0};
    const std::size_t basis_count = harmonics.size();
    if (basis_count > std::numeric_limits<std::size_t>::max() / kShapeSamples ||
        basis_count * kShapeSamples > std::numeric_limits<std::size_t>::max() / kShapeComponents)
        return false;
    out.assign(kShapeSamples * basis_count * kShapeComponents, 0.0);
    for (std::size_t i = 0; i < kShapeSamples; ++i) {
        const double x = static_cast<double>(i) * static_cast<double>(n - 1) /
                         static_cast<double>(kShapeSamples - 1);
        const auto lo = static_cast<std::size_t>(x);
        const auto hi = std::min(lo + 1, n - 1);
        const double t = x - static_cast<double>(lo);
        for (std::size_t b = 0; b < basis_count; ++b) {
            const std::size_t dst = (i * basis_count + b) * kShapeComponents;
            if (has_fourier) {
                const std::size_t a = lo * basis_count + b;
                const std::size_t c = hi * basis_count + b;
                out[dst] = (1.0 - t) * mode.outer_surface_fourier_radial[a] +
                           t * mode.outer_surface_fourier_radial[c];
                out[dst + 1u] = (1.0 - t) * mode.outer_surface_fourier_azimuthal[a] +
                                t * mode.outer_surface_fourier_azimuthal[c];
                out[dst + 2u] = (1.0 - t) * mode.outer_surface_fourier_axial[a] +
                                t * mode.outer_surface_fourier_axial[c];
            } else {
                out[dst] = (1.0 - t) * mode.outer_surface_radial[lo] +
                           t * mode.outer_surface_radial[hi];
                out[dst + 1u] = 0.0;
                out[dst + 2u] = (1.0 - t) * mode.outer_surface_axial[lo] +
                                t * mode.outer_surface_axial[hi];
            }
        }
    }
    double norm2 = 0.0;
    for (double v : out) norm2 += v * v;
    if (!(norm2 > 1e-28) || !std::isfinite(norm2)) return false;
    const double inv = 1.0 / std::sqrt(norm2);
    for (double& v : out) v *= inv;
    return true;
}

std::vector<int> commonHarmonics(const std::vector<int>& a, const std::vector<int>& b,
                                 const std::vector<int>& c, const std::vector<int>& d) {
    std::vector<int> result;
    for (int h : a) {
        if (std::find(b.begin(), b.end(), h) == b.end()) continue;
        if (std::find(c.begin(), c.end(), h) == c.end()) continue;
        if (std::find(d.begin(), d.end(), h) == d.end()) continue;
        if (std::find(result.begin(), result.end(), h) == result.end()) result.push_back(h);
    }
    return result;
}

bool alignShape(const std::vector<double>& source, const std::vector<int>& source_harmonics,
                const std::vector<int>& target_harmonics, std::vector<double>& aligned) {
    if (source_harmonics.empty() || target_harmonics.empty()) return false;
    if (source_harmonics.size() > std::numeric_limits<std::size_t>::max() / kShapeSamples ||
        source_harmonics.size() * kShapeSamples > std::numeric_limits<std::size_t>::max() / kShapeComponents ||
        target_harmonics.size() > std::numeric_limits<std::size_t>::max() / kShapeSamples ||
        target_harmonics.size() * kShapeSamples > std::numeric_limits<std::size_t>::max() / kShapeComponents) return false;
    const std::size_t source_size = kShapeSamples * source_harmonics.size() * kShapeComponents;
    const std::size_t target_size = kShapeSamples * target_harmonics.size() * kShapeComponents;
    if (source.size() != source_size) return false;
    aligned.assign(target_size, 0.0);
    for (std::size_t b = 0; b < target_harmonics.size(); ++b) {
        const auto it = std::find(source_harmonics.begin(), source_harmonics.end(), target_harmonics[b]);
        if (it == source_harmonics.end()) continue; // zero-pad a newly introduced Fourier component
        const std::size_t src_b = static_cast<std::size_t>(it - source_harmonics.begin());
        for (std::size_t i = 0; i < kShapeSamples; ++i) {
            const std::size_t src = (i * source_harmonics.size() + src_b) * kShapeComponents;
            const std::size_t dst = (i * target_harmonics.size() + b) * kShapeComponents;
            for (std::size_t c = 0; c < kShapeComponents; ++c)
                aligned[dst + c] = source[src + c];
        }
    }
    double norm2 = 0.0;
    for (double v : aligned) norm2 += v * v;
    if (!(norm2 > 1e-28)) return false;
    const double inv = 1.0 / std::sqrt(norm2);
    for (double& v : aligned) v *= inv;
    return true;
}

bool orthonormalPair(const std::vector<double>& a, const std::vector<double>& b,
                    std::vector<double>& q0, std::vector<double>& q1) {
    if (a.size() != b.size() || a.empty()) return false;
    q0 = a;
    q1 = b;
    double n0 = 0.0;
    for (double v : q0) n0 += v * v;
    if (!(n0 > 1e-20)) return false;
    const double inv0 = 1.0 / std::sqrt(n0);
    for (double& v : q0) v *= inv0;
    double projection = 0.0;
    for (std::size_t i = 0; i < q0.size(); ++i) projection += q0[i] * q1[i];
    double n1 = 0.0;
    for (std::size_t i = 0; i < q1.size(); ++i) {
        q1[i] -= projection * q0[i];
        n1 += q1[i] * q1[i];
    }
    if (!(n1 > 1e-16)) return false;
    const double inv1 = 1.0 / std::sqrt(n1);
    for (double& v : q1) v *= inv1;
    return true;
}

double subspaceSimilarity(const std::vector<double>& a0, const std::vector<double>& a1,
                          const std::vector<double>& b0, const std::vector<double>& b1) {
    std::vector<double> qa0, qa1, qb0, qb1;
    if (!orthonormalPair(a0, a1, qa0, qa1) ||
        !orthonormalPair(b0, b1, qb0, qb1)) return 0.0;
    double overlap_frobenius2 = 0.0;
    const std::vector<double>* left[2] = {&qa0, &qa1};
    const std::vector<double>* right[2] = {&qb0, &qb1};
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        double dot = 0.0;
        for (std::size_t k = 0; k < left[i]->size(); ++k) dot += (*left[i])[k] * (*right[j])[k];
        overlap_frobenius2 += dot * dot;
    }
    return std::min(1.0, 0.5 * overlap_frobenius2);
}

double shapeSimilarity(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) return 0.0;
    double dot = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) dot += a[i] * b[i];
    return std::min(1.0, std::abs(dot));
}

// Full-volume direct-FEM comparison using a row-sum lumped mass metric.
// Only compares a common physical mesh and common Fourier components; when
// those preconditions are absent the caller deliberately falls back to the
// boundary-trace diagnostic instead of pretending unlike fields are aligned.
bool massWeightedFieldSimilarity(const ModalMode& a, const ModalMode& b, double& similarity) {
    similarity = 0.0;
    const std::size_t na = a.fem_node_r.size();
    const std::size_t nb = b.fem_node_r.size();
    if (na < 2 || na != nb || a.fem_node_z.size() != na || b.fem_node_z.size() != nb ||
        a.fem_harmonics.empty() || b.fem_harmonics.empty()) return false;
    const std::size_t ba = a.fem_harmonics.size(), bb = b.fem_harmonics.size();
    if (na > std::numeric_limits<std::size_t>::max() / ba ||
        nb > std::numeric_limits<std::size_t>::max() / bb) return false;
    const std::size_t a_nodes_basis = na * ba;
    const std::size_t b_nodes_basis = nb * bb;
    if (a_nodes_basis > std::numeric_limits<std::size_t>::max() / 3u ||
        b_nodes_basis > std::numeric_limits<std::size_t>::max() / 3u) return false;
    const std::size_t a_size = a_nodes_basis * 3u;
    const std::size_t b_size = b_nodes_basis * 3u;
    if (a.fem_displacement.size() != a_size || b.fem_displacement.size() != b_size ||
        a.fem_mass_lumped_weights.size() != a_size ||
        b.fem_mass_lumped_weights.size() != b_size) return false;
    for (std::size_t n = 0; n < na; ++n) {
        const double scale = std::max({1.0, std::abs(a.fem_node_r[n]), std::abs(a.fem_node_z[n])});
        if (std::abs(a.fem_node_r[n] - b.fem_node_r[n]) > 1e-9 * scale ||
            std::abs(a.fem_node_z[n] - b.fem_node_z[n]) > 1e-9 * scale) return false;
    }
    double dot = 0.0, norm_a = 0.0, norm_b = 0.0;
    std::size_t common = 0;
    for (std::size_t ha = 0; ha < ba; ++ha) {
        const auto it = std::find(b.fem_harmonics.begin(), b.fem_harmonics.end(), a.fem_harmonics[ha]);
        if (it == b.fem_harmonics.end()) continue;
        const std::size_t hb = static_cast<std::size_t>(it - b.fem_harmonics.begin());
        ++common;
        for (std::size_t n = 0; n < na; ++n) for (std::size_t c = 0; c < 3; ++c) {
            // Direct cylindrical FEM stores [node][component][Fourier basis].
            const std::size_t ia = (n * 3u + c) * ba + ha;
            const std::size_t ib = (n * 3u + c) * bb + hb;
            const double wa = a.fem_mass_lumped_weights[ia];
            const double wb = b.fem_mass_lumped_weights[ib];
            const double w = std::sqrt(std::max(0.0, wa) * std::max(0.0, wb));
            const double va = a.fem_displacement[ia], vb = b.fem_displacement[ib];
            dot += w * va * vb;
            norm_a += w * va * va;
            norm_b += w * vb * vb;
        }
    }
    if (common == 0 || !(norm_a > 1e-28) || !(norm_b > 1e-28)) return false;
    similarity = std::min(1.0, std::abs(dot) / std::sqrt(norm_a * norm_b));
    return std::isfinite(similarity);
}

} // namespace

ModalStateTracker::ModalStateTracker(const ModalTrackerConfig& config) : config_(config) {}

void ModalStateTracker::reset() {
    tracks_ = {};
    last_parameter_ = 0.0;
    initialized_ = false;
}

bool ModalStateTracker::initialize(const ModalData& modes, std::size_t first_index,
                                   std::size_t second_index, double parameter_value) {
    if (first_index >= modes.modes.size() || second_index >= modes.modes.size() ||
        first_index == second_index || !std::isfinite(parameter_value)) return false;
    Track a, b;
    const auto& ma = modes.modes[first_index];
    const auto& mb = modes.modes[second_index];
    if (!(ma.frequency_hz > 0.0) || !(mb.frequency_hz > 0.0) ||
        !makeShape(ma, a.shape, a.harmonics) || !makeShape(mb, b.shape, b.harmonics)) return false;
    a.frequency = ma.frequency_hz;
    b.frequency = mb.frequency_hz;
    a.full_mode = ma;
    b.full_mode = mb;
    a.p00 = b.p00 = std::max(1.0, config_.measurement_frequency_variance);
    a.p11 = b.p11 = 100.0;
    tracks_[0] = std::move(a);
    tracks_[1] = std::move(b);
    last_parameter_ = parameter_value;
    initialized_ = true;
    return true;
}

ModalStateTrackerReport ModalStateTracker::update(const ModalData& candidates,
                                                   double parameter_value) {
    ModalStateTrackerReport report;
    report.initialized = initialized_;
    if (!initialized_ || candidates.modes.size() < 2 || !std::isfinite(parameter_value)) return report;

    const double dt = parameter_value - last_parameter_;
    if (std::abs(dt) < 1e-15) return report;

    struct Predicted { double f, slope, p00, p01, p11, sigma; } pred[2];
    for (int i = 0; i < 2; ++i) {
        const auto& t = tracks_[static_cast<std::size_t>(i)];
        auto& p = pred[i];
        p.f = t.frequency + dt * t.slope;
        p.slope = t.slope;
        const double qf = config_.process_frequency_variance * std::abs(dt);
        const double qs = config_.process_slope_variance * std::abs(dt);
        p.p00 = t.p00 + 2.0 * dt * t.p01 + dt * dt * t.p11 + qf;
        p.p01 = t.p01 + dt * t.p11;
        p.p11 = t.p11 + qs;
        p.sigma = std::sqrt(std::max(1e-12, p.p00 + config_.measurement_frequency_variance));
    }

    struct Hypothesis { std::size_t a, b; double score; double sim_a, sim_b, pair_similarity; };
    std::vector<std::vector<double>> shapes(candidates.modes.size());
    std::vector<std::vector<int>> shape_harmonics(candidates.modes.size());
    std::vector<bool> valid(candidates.modes.size(), false);
    for (std::size_t j = 0; j < candidates.modes.size(); ++j) {
        const auto& m = candidates.modes[j];
        valid[j] = m.frequency_hz > 0.0 && std::isfinite(m.frequency_hz) &&
                   makeShape(m, shapes[j], shape_harmonics[j]);
    }
    std::vector<Hypothesis> hypotheses;
    for (std::size_t a = 0; a < candidates.modes.size(); ++a) {
        if (!valid[a]) continue;
        const double da = std::abs(std::log(candidates.modes[a].frequency_hz / pred[0].f));
        if (da > config_.maximum_log_frequency_error) continue;
        for (std::size_t b = 0; b < candidates.modes.size(); ++b) {
            if (a == b || !valid[b]) continue;
            const double db = std::abs(std::log(candidates.modes[b].frequency_hz / pred[1].f));
            if (db > config_.maximum_log_frequency_error) continue;
            // Compare only the Fourier components present in both the tracked
            // family and the candidate family.  A mesh/Fourier enrichment must
            // not look like a physical branch change merely because the new
            // basis contains additional (initially untracked) components.
            const auto comparison_harmonics = commonHarmonics(
                tracks_[0].harmonics, tracks_[1].harmonics,
                shape_harmonics[a], shape_harmonics[b]);
            if (comparison_harmonics.empty()) continue;
            std::vector<double> ta, tb, ca, cb;
            if (!alignShape(tracks_[0].shape, tracks_[0].harmonics, comparison_harmonics, ta) ||
                !alignShape(tracks_[1].shape, tracks_[1].harmonics, comparison_harmonics, tb) ||
                !alignShape(shapes[a], shape_harmonics[a], comparison_harmonics, ca) ||
                !alignShape(shapes[b], shape_harmonics[b], comparison_harmonics, cb)) continue;
            ++report.frequency_compatible_assignments;
            const double pair_similarity = subspaceSimilarity(ta, tb, ca, cb);
            if (pair_similarity < config_.minimum_subspace_similarity) {
                ++report.subspace_rejected_assignments;
                report.best_rejected_subspace_similarity = std::max(
                    report.best_rejected_subspace_similarity, pair_similarity);
                continue;
            }
            const double fa = std::exp(-0.5 * da * da * pred[0].f * pred[0].f /
                                       (pred[0].sigma * pred[0].sigma));
            const double fb = std::exp(-0.5 * db * db * pred[1].f * pred[1].f /
                                       (pred[1].sigma * pred[1].sigma));
            double sa = shapeSimilarity(ta, ca);
            double sb = shapeSimilarity(tb, cb);
            double full_sa = 0.0, full_sb = 0.0;
            if (massWeightedFieldSimilarity(tracks_[0].full_mode, candidates.modes[a], full_sa)) sa = full_sa;
            if (massWeightedFieldSimilarity(tracks_[1].full_mode, candidates.modes[b], full_sb)) sb = full_sb;
            // The pair subspace identifies the family, but is invariant to a
            // swap of its two basis vectors. Individual shape continuity must
            // therefore contribute explicitly to distinguish the two branches.
            // Frequency prediction remains a weaker tie-breaker near crossings.
            const double score = 0.35 * pair_similarity +
                                 0.25 * (sa + sb) +
                                 0.075 * (fa + fb);
            hypotheses.push_back({a, b, score, sa, sb, pair_similarity});
        }
    }
    if (hypotheses.empty()) {
        report.no_compatible_assignment = true;
        return report;
    }
    std::sort(hypotheses.begin(), hypotheses.end(), [](const auto& x, const auto& y) {
        return x.score > y.score;
    });
    const auto& best = hypotheses.front();
    report.best_assignment_score = best.score;
    report.pair_subspace_similarity = best.pair_similarity;
    report.assignment_margin = hypotheses.size() > 1 ? best.score - hypotheses[1].score : best.score;
    report.ambiguous = hypotheses.size() > 1 &&
                       report.assignment_margin < config_.minimum_assignment_margin;
    for (int i = 0; i < 2; ++i) {
        const std::size_t index = i == 0 ? best.a : best.b;
        auto& e = report.tracks[static_cast<std::size_t>(i)];
        e.valid = true;
        e.candidate_index = index;
        e.predicted_frequency_hz = pred[i].f;
        e.measured_frequency_hz = candidates.modes[index].frequency_hz;
        e.frequency_sigma_hz = pred[i].sigma;
        e.shape_similarity = i == 0 ? best.sim_a : best.sim_b;
        e.assignment_score = best.score;
    }
    if (report.ambiguous) return report;

    for (int i = 0; i < 2; ++i) {
        auto& t = tracks_[static_cast<std::size_t>(i)];
        const auto& p = pred[i];
        const auto& e = report.tracks[static_cast<std::size_t>(i)];
        const double r = std::max(1e-12, config_.measurement_frequency_variance);
        const double innovation = e.measured_frequency_hz - p.f;
        const double denom = p.p00 + r;
        const double k0 = p.p00 / denom;
        const double k1 = p.p01 / denom;
        t.frequency = p.f + k0 * innovation;
        t.slope = p.slope + k1 * innovation;
        t.p00 = std::max(1e-12, (1.0 - k0) * p.p00);
        t.p01 = (1.0 - k0) * p.p01;
        t.p11 = std::max(1e-12, p.p11 - k1 * p.p01);
        t.shape = shapes[e.candidate_index];
        t.harmonics = shape_harmonics[e.candidate_index];
        t.full_mode = candidates.modes[e.candidate_index];
    }
    last_parameter_ = parameter_value;
    report.updated = true;
    return report;
}

} // namespace BellOptimizer
