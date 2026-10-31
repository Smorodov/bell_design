#include "modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace BellOptimizer {
namespace {

bool hasShape(const ModalMode& mode) {
    return !mode.outer_surface_radial.empty() &&
           mode.outer_surface_radial.size() == mode.outer_surface_axial.size();
}

double shapeCosine(const ModalMode& a, const ModalMode& b) {
    if (!hasShape(a) || !hasShape(b) ||
        a.outer_surface_radial.size() != b.outer_surface_radial.size()) {
        return -1.0;
    }
    double dot = 0.0, na = 0.0, nb = 0.0;
    for (std::size_t i = 0; i < a.outer_surface_radial.size(); ++i) {
        const double ar = a.outer_surface_radial[i];
        const double az = a.outer_surface_axial[i];
        const double br = b.outer_surface_radial[i];
        const double bz = b.outer_surface_axial[i];
        dot += ar * br + az * bz;
        na += ar * ar + az * az;
        nb += br * br + bz * bz;
    }
    if (!(na > 1e-30) || !(nb > 1e-30)) return -1.0;
    return std::abs(dot) / std::sqrt(na * nb);
}

// Small rectangular assignment solved by dynamic programming over subsets.
// The number of tracked modes is intentionally capped by the available modal
// result (normally 10), so O(N 2^N) remains negligible compared with FEM.
std::vector<int> maximumAssignment(const std::vector<std::vector<double>>& score) {
    const std::size_t n = score.size();
    if (n == 0) return {};
    const std::size_t m = score.front().size();
    if (m == 0) return std::vector<int>(n, -1);

    const std::size_t rows = std::min(n, m);
    const std::size_t cols = std::max(n, m);

    if (cols > 20) {
        std::vector<int> result(n, -1);
        std::vector<bool> used(m, false);
        for (std::size_t i = 0; i < n; ++i) {
            double best = -1e9;
            int best_j = -1;
            for (std::size_t j = 0; j < m; ++j) {
                if (!used[j] && score[i][j] > best) {
                    best = score[i][j];
                    best_j = static_cast<int>(j);
                }
            }
            if (best_j >= 0) { result[i] = best_j; used[best_j] = true; }
        }
        return result;
    }

    const std::size_t states = std::size_t(1) << cols;
    std::vector<double> dp(states, -std::numeric_limits<double>::infinity());
    std::vector<std::vector<int>> parent(rows, std::vector<int>(states, -1));
    dp[0] = 0.0;

    for (std::size_t i = 0; i < rows; ++i) {
        std::vector<double> next(states, -std::numeric_limits<double>::infinity());
        for (std::size_t mask = 0; mask < states; ++mask) {
            if (!std::isfinite(dp[mask])) continue;
            for (std::size_t j = 0; j < cols; ++j) {
                if (mask & (std::size_t(1) << j)) continue;
                const std::size_t nm = mask | (std::size_t(1) << j);
                const double value = dp[mask] + (n <= m ? score[i][j] : score[j][i]);
                if (value > next[nm]) {
                    next[nm] = value;
                    parent[i][nm] = static_cast<int>(j);
                }
            }
        }
        dp.swap(next);
    }

    std::size_t mask = 0;
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t s = 0; s < states; ++s) {
        if (std::isfinite(dp[s]) && dp[s] > best) { best = dp[s]; mask = s; }
    }

    std::vector<int> compact(rows, -1);
    for (std::size_t i = rows; i-- > 0;) {
        const int j = parent[i][mask];
        compact[i] = j;
        if (j >= 0) mask &= ~(std::size_t(1) << static_cast<std::size_t>(j));
    }

    if (n <= m) return compact;

    std::vector<int> result(n, -1);
    for (std::size_t ref = 0; ref < compact.size(); ++ref) {
        const int candidate = compact[ref];
        if (candidate >= 0 && static_cast<std::size_t>(candidate) < result.size())
            result[static_cast<std::size_t>(candidate)] = static_cast<int>(ref);
    }
    return result;
}
} // namespace

ModalTrackingReport trackModalFamilies(const ModalData& reference,
                                       ModalData& candidate,
                                       double ambiguity_similarity,
                                       double near_degenerate_frequency_ratio) {
    ModalTrackingReport report;
    if (reference.modes.empty() || candidate.modes.empty()) return report;

    const std::size_t n = reference.modes.size();
    const std::size_t m = candidate.modes.size();
    const std::size_t count = std::min<std::size_t>({n, m, 10});
    if (count < 2) return report;

    for (std::size_t i = 0; i < count; ++i) {
        if (!hasShape(reference.modes[i])) return report;
    }
    for (std::size_t j = 0; j < count; ++j) {
        if (!hasShape(candidate.modes[j])) return report;
    }

    std::vector<std::vector<double>> score(count, std::vector<double>(count, -1.0));
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            const double shape = shapeCosine(reference.modes[i], candidate.modes[j]);
            if (shape < 0.0) continue;
            // Shape is the identity criterion. Frequency is only a weak tie-breaker,
            // so a genuine crossing does not silently relabel the modal family.
            double frequency_term = 0.0;
            if (reference.modes[i].frequency_hz > 0.0 && candidate.modes[j].frequency_hz > 0.0) {
                const double ratio = candidate.modes[j].frequency_hz / reference.modes[i].frequency_hz;
                if (ratio > 0.0 && std::isfinite(ratio)) {
                    frequency_term = std::exp(-std::abs(std::log(ratio)));
                }
            }
            score[i][j] = 0.90 * shape + 0.10 * frequency_term;
        }
    }

    // Here rows are reference modes and columns are candidate modes. For the
    // square case the assignment directly gives candidate index per reference.
    const auto assignment = maximumAssignment(score);
    if (assignment.size() < count) return report;

    std::vector<ModalMode> tracked;
    tracked.reserve(m);
    std::vector<bool> used(m, false);
    std::size_t matched = 0;
    double min_similarity = 1.0;
    double sum_similarity = 0.0;

    for (std::size_t ref = 0; ref < count; ++ref) {
        const int cj = assignment[ref];
        if (cj < 0 || static_cast<std::size_t>(cj) >= m || used[static_cast<std::size_t>(cj)]) continue;
        const auto& source = candidate.modes[static_cast<std::size_t>(cj)];
        const double similarity = shapeCosine(reference.modes[ref], source);
        if (!(similarity >= 0.0)) continue;
        tracked.push_back(source);
        used[static_cast<std::size_t>(cj)] = true;
        ++matched;
        min_similarity = std::min(min_similarity, similarity);
        sum_similarity += similarity;

        const double rf = reference.modes[ref].frequency_hz;
        const double cf = source.frequency_hz;
        if (rf > 0.0 && cf > 0.0) {
            const double frequency_ratio = std::abs(cf / rf - 1.0);
            if (frequency_ratio < near_degenerate_frequency_ratio && similarity < 0.75)
                report.ambiguous = true;
        }
    }

    if (matched < count) return report;

    // Keep all untracked higher modes after the tracked modal families, in their
    // original solver order. This preserves deterministic behaviour without
    // re-sorting away the family identity we just established.
    for (std::size_t j = 0; j < m; ++j) {
        if (!used[j]) tracked.push_back(candidate.modes[j]);
    }
    candidate.modes = std::move(tracked);

    report.used_shape_tracking = true;
    report.matched_modes = matched;
    report.minimum_similarity = min_similarity;
    report.mean_similarity = sum_similarity / static_cast<double>(matched);
    if (report.minimum_similarity < ambiguity_similarity) report.ambiguous = true;
    return report;
}

} // namespace BellOptimizer
