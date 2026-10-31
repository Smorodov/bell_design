#include "fem/sparse_eigensolver.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    using BellOptimizer::detail::SparseMatrix;
    using BellOptimizer::detail::SparseTriplet;
    using BellOptimizer::detail::solve_lowest_sparse_modes;
    const SparseMatrix k = SparseMatrix::from_triplets(5, {
        {0,0,1.0}, {1,1,4.0}, {2,2,9.0}, {3,3,16.0}, {4,4,25.0}
    });
    const SparseMatrix m = SparseMatrix::from_triplets(5, {
        {0,0,1.0}, {1,1,1.0}, {2,2,1.0}, {3,3,1.0}, {4,4,1.0}
    });
    const auto result = solve_lowest_sparse_modes(k, m, 2, 100, 1e-8);
    if (!result.converged || result.eigenvalues.size() != 2 || result.vectors.size() != 2) return 1;
    if (std::abs(result.eigenvalues[0] - 1.0) > 1e-6 ||
        std::abs(result.eigenvalues[1] - 4.0) > 1e-5) return 2;
    for (std::size_t mode = 0; mode < result.vectors.size(); ++mode) {
        const auto& v = result.vectors[mode];
        std::vector<double> residual(5, 0.0);
        for (int i = 0; i < 5; ++i) residual[static_cast<std::size_t>(i)] =
            (i + 1) * (i + 1) * v[static_cast<std::size_t>(i)] -
            result.eigenvalues[mode] * v[static_cast<std::size_t>(i)];
        double norm = 0.0; for (double x : residual) norm += x*x;
        if (std::sqrt(norm) > 1e-5) return 3;
    }
    std::cout << "sparse_eigensolver_contract=PASS\n"
              << "iterations=" << result.iterations << '\n'
              << "lambda0=" << result.eigenvalues[0] << '\n'
              << "lambda1=" << result.eigenvalues[1] << '\n';
    return 0;
}
