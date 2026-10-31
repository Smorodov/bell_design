#include "fem/sparse_eigensolver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;
using BellOptimizer::detail::solve_lowest_sparse_modes;

static SparseMatrix diagonal(std::initializer_list<double> values) {
    std::vector<SparseTriplet> entries;
    int i = 0;
    for (double value : values) entries.push_back({i, i++, value});
    return SparseMatrix::from_triplets(i, std::move(entries));
}

int main() {
    const auto k = diagonal({1.0, 4.0, 9.0});
    const auto m = diagonal({1.0, 1.0, 1.0});

    const auto good = solve_lowest_sparse_modes(k, m, 2, 100, 1e-8, true);
    if (!good.converged || good.eigenvalues.size() != 2 ||
        !std::isfinite(good.final_relative_residual) ||
        good.final_relative_residual >= 1e-8) return 1;

    auto nonsymmetric = SparseMatrix::from_triplets(3, {
        {0,0,2.0}, {0,1,-1.0}, {1,0,-1.0001}, {1,1,2.0}, {2,2,2.0}
    });
    if (solve_lowest_sparse_modes(nonsymmetric, m, 1, 20, 1e-8).converged) return 2;

    // A missing upper-triangular counterpart must also be rejected; checking
    // only the upper triangle would incorrectly accept this matrix as symmetric.
    auto one_sided = SparseMatrix::from_triplets(3, {
        {0,0,2.0}, {1,0,-1.0}, {1,1,2.0}, {2,2,2.0}
    });
    if (solve_lowest_sparse_modes(one_sided, m, 1, 20, 1e-8).converged) return 3;

    auto nonfinite = SparseMatrix::from_triplets(3, {
        {0,0,std::numeric_limits<double>::quiet_NaN()}, {1,1,2.0}, {2,2,2.0}
    });
    if (solve_lowest_sparse_modes(nonfinite, m, 1, 20, 1e-8).converged) return 4;

    const auto bad_mass = diagonal({1.0, 0.0, 1.0});
    if (solve_lowest_sparse_modes(k, bad_mass, 1, 20, 1e-8).converged) return 5;
    if (solve_lowest_sparse_modes(k, m, 1, 0, 1e-8).converged) return 6;
    if (solve_lowest_sparse_modes(k, m, 1, 20, 0.0).converged) return 7;

    if (std::abs(good.eigenvalues[0] - 1.0) > 1e-7 ||
        std::abs(good.eigenvalues[1] - 4.0) > 1e-6 ||
        good.max_mass_orthogonality_error >= 1e-8 ||
        good.max_eigenpair_relative_residual >= 1e-8 ||
        !std::isfinite(good.regularization_shift) || good.regularization_shift <= 0.0) return 8;

    std::cout << "sparse_solver_production_contract=PASS\n"
              << "residual=" << good.final_relative_residual << '\n'
              << "pcg_iterations=" << good.pcg_iterations << '\n'
              << "mass_orthogonality_error=" << good.max_mass_orthogonality_error << '\n'
              << "eigenpair_residual=" << good.max_eigenpair_relative_residual << '\n'
              << "regularization_shift=" << good.regularization_shift << '\n';
    return 0;
}
