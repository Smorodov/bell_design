#include "fem/sparse_eigensolver.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;
using BellOptimizer::detail::solve_lowest_sparse_modes;

static SparseMatrix make_stiffness(int n) {
    std::vector<SparseTriplet> entries;
    entries.reserve(static_cast<std::size_t>(3 * n - 2));
    for (int i = 0; i < n; ++i) {
        entries.push_back({i, i, 2.0});
        if (i > 0) {
            entries.push_back({i, i - 1, -1.0});
            entries.push_back({i - 1, i, -1.0});
        }
    }
    return SparseMatrix::from_triplets(n, std::move(entries));
}

static SparseMatrix make_mass(int n) {
    std::vector<SparseTriplet> entries;
    entries.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) entries.push_back({i, i, 1.0});
    return SparseMatrix::from_triplets(n, std::move(entries));
}

int main() {
    for (const int n : {1000, 2000, 5000}) {
        const auto k = make_stiffness(n);
        const auto m = make_mass(n);
        const auto start = std::chrono::steady_clock::now();
        const auto result = solve_lowest_sparse_modes(k, m, 5, 100, 1e-8, true);
        const double elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        if (!result.converged || result.eigenvalues.size() != 5 ||
            result.vectors.size() != 5 ||
            !std::isfinite(result.final_relative_residual) ||
            result.final_relative_residual >= 1e-8) return 1;

        const double exact_lambda0 = 2.0 - 2.0 * std::cos(std::acos(-1.0) / static_cast<double>(n + 1));
        const double relative_error = std::abs(result.eigenvalues[0] - exact_lambda0) /
            std::max(std::abs(exact_lambda0), 1e-300);
        if (relative_error > 1e-5 || result.pcg_iterations > 2000) return 2;

        std::cout << "scale_case=PASS n=" << n
                  << " nnz=" << k.values.size()
                  << " ic0_nnz=" << result.ic0_nonzeros
                  << " outer_iterations=" << result.iterations
                  << " pcg_iterations=" << result.pcg_iterations
                  << " residual=" << result.final_relative_residual
                  << " lambda0_rel_error=" << relative_error
                  << " elapsed_ms=" << elapsed_ms << '\n';
    }
    std::cout << "sparse_solver_scaling=PASS\n";
    return 0;
}
