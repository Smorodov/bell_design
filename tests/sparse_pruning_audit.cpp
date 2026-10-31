#include "fem/sparse_matrix.hpp"
#include "fem/sparse_eigensolver.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;
using BellOptimizer::detail::solve_lowest_sparse_modes;
using BellOptimizer::detail::prune_symmetric_offdiagonals;

static double max_relative_eigenvalue_change(const std::vector<double>& a,
                                             const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        e = std::max(e, std::abs(a[i] - b[i]) /
            std::max(std::abs(a[i]), std::numeric_limits<double>::min()));
    return e;
}

int main() {
    // SPD chain with deliberately weak long-range couplings. The weak links
    // are candidates for removal; nearest-neighbour links are not.
    constexpr int n = 40;
    std::vector<SparseTriplet> k_entries;
    std::vector<SparseTriplet> m_entries;
    for (int i = 0; i < n; ++i) {
        k_entries.push_back({i, i, 2.0 + (i == 0 || i == n - 1 ? -1.0 : 0.0)});
        m_entries.push_back({i, i, 1.0});
        if (i + 1 < n) {
            k_entries.push_back({i, i + 1, -1.0});
            k_entries.push_back({i + 1, i, -1.0});
        }
        if (i + 3 < n) {
            k_entries.push_back({i, i + 3, 2.0e-7});
            k_entries.push_back({i + 3, i, 2.0e-7});
        }
    }
    const auto k = SparseMatrix::from_triplets(n, std::move(k_entries));
    const auto m = SparseMatrix::from_triplets(n, std::move(m_entries));

    const auto baseline = solve_lowest_sparse_modes(k, m, 5, 160, 1e-8);
    if (!baseline.converged || baseline.eigenvalues.size() < 5) {
        std::cerr << "baseline sparse solve did not converge\n";
        return 1;
    }

    const double thresholds[] = {1e-7, 1e-6, 1e-5, 1e-4};
    std::size_t previous_removed = 0;
    for (double threshold : thresholds) {
        const auto report = prune_symmetric_offdiagonals(k, threshold);
        const auto& pruned = report.matrix;
        if (report.removed_pairs < previous_removed) {
            std::cerr << "pruning is not monotone\n";
            return 2;
        }
        previous_removed = report.removed_pairs;
        if (report.removed_pairs == 0 || pruned.values.size() >= k.values.size()) {
            std::cerr << "expected weak symmetric pairs to be removed at threshold " << threshold << "\n";
            return 3;
        }
        const double norm_bound = report.original_frobenius_norm > 0.0
            ? report.removed_frobenius_norm / report.original_frobenius_norm : 0.0;
        if (!(norm_bound < 1e-5)) {
            std::cerr << "removed norm unexpectedly large: " << norm_bound << "\n";
            return 4;
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (std::abs(pruned.get(i,j) - pruned.get(j,i)) > 1e-14) {
                    std::cerr << "pruned matrix lost symmetry\n";
                    return 5;
                }
        const auto result = solve_lowest_sparse_modes(pruned, m, 5, 160, 1e-8);
        if (!result.converged || result.eigenvalues.size() < 5) {
            std::cerr << "pruned solve failed at threshold " << threshold << "\n";
            return 6;
        }
        const double spectral_change = max_relative_eigenvalue_change(
            baseline.eigenvalues, result.eigenvalues);
        std::cout << "threshold=" << threshold
                  << " removed_pairs=" << report.removed_pairs
                  << " nnz=" << pruned.values.size()
                  << " removed_frobenius_ratio=" << norm_bound
                  << " max_relative_eigenvalue_change=" << spectral_change << '\n';
        if (!std::isfinite(spectral_change)) {
            std::cerr << "non-finite spectral change at threshold " << threshold << "\n";
            return 7;
        }
        // The first mode is deliberately near-null. This is an important
        // negative result: a tiny Frobenius perturbation can still cause a
        // large relative shift of a low eigenvalue. At the first threshold
        // the audit therefore requires the large change to remain visible.
        if (threshold == 1e-7 && !(spectral_change > 0.1)) {
            std::cerr << "near-null sensitivity was not detected\n";
            return 8;
        }
    }

    // A non-symmetric pair must never be silently removed.
    const auto nonsym = SparseMatrix::from_triplets(2, {
        {0,0,2.0}, {1,1,2.0}, {0,1,1e-9}, {1,0,2e-9}
    });
    const auto nonsym_report = prune_symmetric_offdiagonals(nonsym, 1e-4);
    if (nonsym_report.removed_pairs != 0 || nonsym_report.matrix.values.size() != 4u) {
        std::cerr << "non-symmetric weak pair was incorrectly pruned\n";
        return 9;
    }

    std::cout << "sparse_pruning_audit=PASS\n";
    return 0;
}
