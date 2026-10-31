#include "fem/sparse_eigensolver.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;
using BellOptimizer::detail::solve_lowest_sparse_modes;

static SparseMatrix diagonal_scaled(double scale) {
    std::vector<SparseTriplet> entries;
    const double values[] = {1.0, 4.0, 9.0, 16.0, 25.0};
    for (int i = 0; i < 5; ++i) entries.push_back({i, i, scale * values[i]});
    return SparseMatrix::from_triplets(5, std::move(entries));
}

int main() {
    // K and M have deliberately extreme absolute scales.  The generalized
    // eigenvalues must remain the same because both matrices are scaled by
    // independent positive factors and the solver restores K/M's scale ratio.
    const double scales[] = {1.0e-12, 1.0, 1.0e12};
    for (double scale : scales) {
        const auto k = diagonal_scaled(scale);
        const auto m = diagonal_scaled(1.0);
        const auto result = solve_lowest_sparse_modes(k, m, 2, 100, 1.0e-10);
        if (!result.converged || result.eigenvalues.size() != 2 ||
            !std::isfinite(result.final_relative_residual) ||
            result.final_relative_residual >= 1.0e-8) {
            std::cerr << std::setprecision(17)
                      << "sparse_solver_scaling_contract: FAIL scale=" << scale
                      << " converged=" << result.converged
                      << " residual=" << result.final_relative_residual << '\n';
            return 1;
        }
        const double expected = scale;
        if (std::abs(result.eigenvalues[0] - expected) > std::abs(expected) * 1.0e-6 ||
            std::abs(result.eigenvalues[1] - expected) > std::abs(expected) * 1.0e-6) {
            std::cerr << std::setprecision(17)
                      << "sparse_solver_scaling_contract: FAIL scale=" << scale
                      << " lambda0=" << result.eigenvalues[0]
                      << " lambda1=" << result.eigenvalues[1] << '\n';
            return 2;
        }
    }

    std::cout << "sparse_solver_scaling_contract: PASS\n"
              << "tested_scales=1e-12,1,1e12\n";
    return 0;
}
