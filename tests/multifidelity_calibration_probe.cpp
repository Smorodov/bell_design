#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "fem/internal_modal_solver.hpp"
#include "core/modal_tracking.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using namespace BellOptimizer;

namespace {
struct CandidateResult {
    std::string name;
    std::array<double, 5> fast{};
    std::array<double, 5> accurate{};
    double fast_score = 0.0;
    double accurate_score = 0.0;
    double minimum_shape_similarity = 0.0;
    double mean_shape_similarity = 0.0;
    double matched_frequency_error = 0.0;
    bool tracking_ambiguous = false;
};

BellGeometry base_geometry() {
    return BellGeometry{658.3175709034034, 355.9758748202059,
        22.85675347179014, 46.62524667849419,
        {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
}

std::vector<std::pair<std::string, BellGeometry>> candidates() {
    const auto base = base_geometry();
    std::vector<std::pair<std::string, BellGeometry>> out;
    out.reserve(30);
    // Deterministic low-amplitude coverage of the three interior profile knots.
    // The perturbations are intentionally modest: this is a local calibration
    // set, not evidence that the entire optimizer search domain is covered.
    constexpr double golden_angle = 2.39996322972865332;
    for (int i = 0; i < 30; ++i) {
        const double p = golden_angle * static_cast<double>(i);
        auto g = base;
        g.r_knots[1] += 0.0040 * std::sin(p);
        g.r_knots[2] += 0.0050 * std::cos(1.31 * p);
        g.r_knots[3] += 0.0060 * std::sin(0.73 * p + 0.4);
        out.emplace_back("candidate_" + std::to_string(i + 1), std::move(g));
    }
    return out;
}

double frequency_shape_score(const std::array<double, 5>& f) {
    if (!(f[1] > 0.0)) return 1.0e9;
    const std::array<double, 4> target{0.5, 1.25, 1.5, 2.0};
    const std::array<int, 4> index{0, 2, 3, 4};
    double score = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const double ratio = f[static_cast<std::size_t>(index[i])] / f[1];
        if (!(ratio > 0.0) || !std::isfinite(ratio)) return 1.0e9;
        const double e = std::log(ratio / target[i]);
        score += e * e;
    }
    return score;
}

bool solve_fast(const BellGeometry& geometry, const BellMaterial& material,
                const BellMeshConfig& mesh_config, std::array<double, 5>& frequencies,
                ModalData* captured = nullptr) {
    BellGenerator generator(material, mesh_config);
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    Excitation::MeridionalInfluence influence;
    ModalProblem problem{generator.getMeshData(), material, ".", "mf_fast", &influence};
    InternalModalSolver solver(true, 8, true);
    ModalData data;
    if (!solver.solve(problem, data) || data.modes.size() < frequencies.size()) return false;
    for (std::size_t i = 0; i < frequencies.size(); ++i)
        frequencies[i] = data.modes[i].frequency_hz;
    if (captured) *captured = std::move(data);
    return true;
}

bool solve_accurate(const BellGeometry& geometry, const BellMaterial& material,
                    const BellMeshConfig& mesh_config, std::array<double, 5>& frequencies,
                    ModalData* captured = nullptr) {
    BellGenerator generator(material, mesh_config);
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh_config);
    Excitation::MeridionalInfluence provider;
    const auto influence = provider.make(model.mesh_si);
    FEM::CylindricalModalEvaluationConfig config;
    // m=0 only: isolate the shared axisymmetric subproblem. This probe does
    // not yet claim equivalence to the complete multi-harmonic acoustic model.
    config.harmonics = {0};
    config.angular_samples = 32;
    config.requested_modes = 8;
    config.max_iterations = 1600;
    config.tolerance = 1.0e-6;
    config.minimum_frequency_hz = 1.0;
    ModalData data;
    if (!FEM::evaluate_cylindrical_bell_modes(model, influence, config, data) ||
        data.modes.size() < frequencies.size()) return false;
    for (std::size_t i = 0; i < frequencies.size(); ++i)
        frequencies[i] = data.modes[i].frequency_hz;
    if (captured) *captured = std::move(data);
    return true;
}

std::vector<double> ranks(const std::vector<double>& values) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return values[a] < values[b]; });
    std::vector<double> out(values.size());
    for (std::size_t i = 0; i < order.size(); ++i) out[order[i]] = static_cast<double>(i + 1);
    return out;
}

double spearman(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.size() < 2) return 0.0;
    const auto ra = ranks(a), rb = ranks(b);
    double ma = 0.0, mb = 0.0;
    for (double x : ra) ma += x;
    for (double x : rb) mb += x;
    ma /= static_cast<double>(ra.size()); mb /= static_cast<double>(rb.size());
    double cov = 0.0, va = 0.0, vb = 0.0;
    for (std::size_t i = 0; i < ra.size(); ++i) {
        const double x = ra[i] - ma, y = rb[i] - mb;
        cov += x * y; va += x * x; vb += y * y;
    }
    return (va > 0.0 && vb > 0.0) ? cov / std::sqrt(va * vb) : 0.0;
}
}

int main() {
    BellMaterial material;
    BellMeshConfig mesh_config; mesh_config.nNodes = 15;
    std::vector<CandidateResult> results;
    for (const auto& item : candidates()) {
        CandidateResult r; r.name = item.first;
        ModalData fast_data, accurate_data;
        if (!solve_fast(item.second, material, mesh_config, r.fast, &fast_data)) {
            std::cerr << "fast evaluator failed: " << r.name << '\n'; return 2;
        }
        if (!solve_accurate(item.second, material, mesh_config, r.accurate, &accurate_data)) {
            std::cerr << "accurate evaluator failed: " << r.name << '\n'; return 3;
        }
        const auto tracking = trackModalFamilies(fast_data, accurate_data, 0.10, 0.01);
        if (!tracking.used_shape_tracking || tracking.matched_modes < 5) {
            std::cerr << "cross-solver shape matching unavailable: " << r.name
                      << " (fast modes=" << fast_data.modes.size()
                      << ", accurate modes=" << accurate_data.modes.size() << ")\n";
            return 7;
        }
        r.minimum_shape_similarity = tracking.minimum_similarity;
        r.mean_shape_similarity = tracking.mean_similarity;
        // A low minimum trace overlap is not a defensible family assignment,
        // even if the legacy tracker can still produce a complete permutation.
        r.tracking_ambiguous = tracking.ambiguous || tracking.minimum_similarity < 0.75;
        for (std::size_t i = 0; i < 5; ++i) {
            r.fast[i] = fast_data.modes[i].frequency_hz;
            r.accurate[i] = accurate_data.modes[i].frequency_hz;
            const double rel = std::abs(r.fast[i] - r.accurate[i]) /
                std::max(std::abs(r.accurate[i]), 1e-12);
            r.matched_frequency_error = std::max(r.matched_frequency_error, rel);
        }
        r.fast_score = frequency_shape_score(r.fast);
        r.accurate_score = frequency_shape_score(r.accurate);
        results.push_back(r);
    }

    std::vector<double> fast_scores, accurate_scores;
    double max_matched_frequency_error = 0.0;
    double min_shape_similarity = 1.0, sum_mean_shape = 0.0;
    std::size_t ambiguous_count = 0;
    std::cout << std::setprecision(9);
    std::cout << "candidate,fast_score,shape_matched_accurate_score,min_shape_similarity,mean_shape_similarity,matched_frequency_error,tracking_ambiguous\n";
    for (const auto& r : results) {
        fast_scores.push_back(r.fast_score); accurate_scores.push_back(r.accurate_score);
        max_matched_frequency_error = std::max(max_matched_frequency_error, r.matched_frequency_error);
        min_shape_similarity = std::min(min_shape_similarity, r.minimum_shape_similarity);
        sum_mean_shape += r.mean_shape_similarity;
        if (r.tracking_ambiguous) ++ambiguous_count;
        std::cout << r.name << ',' << r.fast_score << ',' << r.accurate_score << ','
                  << r.minimum_shape_similarity << ',' << r.mean_shape_similarity << ','
                  << r.matched_frequency_error << ',' << (r.tracking_ambiguous ? 1 : 0) << '\n';
    }
    const double rho = spearman(fast_scores, accurate_scores);
    std::cout << "spearman_shape_matched_objective_rank=" << rho << '\n'
              << "max_matched_first5_frequency_relative_difference=" << max_matched_frequency_error << '\n'
              << "minimum_outer_trace_shape_similarity=" << min_shape_similarity << '\n'
              << "mean_candidate_shape_similarity=" << (sum_mean_shape / static_cast<double>(results.size())) << '\n'
              << "ambiguous_candidate_count=" << ambiguous_count << '\n'
              << "provisionally_usable_candidate_count=" << (results.size() - ambiguous_count) << '\n'
              << "candidate_count=" << results.size() << '\n'
              << "comparison_scope=m0_only_same_15_node_mesh_outer_wall_trace_matching\n";
    // This remains a diagnostic, not a promotion gate. Matching by boundary
    // trace is a useful association hypothesis, not a full-volume proof.
    return std::isfinite(rho) && std::isfinite(max_matched_frequency_error) &&
           std::isfinite(min_shape_similarity) ? 0 : 8;
}
