#include "fem/sparse_eigensolver.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer::detail;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    // A small SPD generalized eigenproblem with deliberately nontrivial graph
    // ordering. The same physical seed must produce a valid continuation solve
    // with and without RCM ordering.
    const std::vector<SparseTriplet> k_entries{
        {0,0,4.0},{0,1,-1.0},{1,0,-1.0},{1,1,4.0},
        {1,2,-1.0},{2,1,-1.0},{2,2,4.0},{2,3,-1.0},
        {3,2,-1.0},{3,3,4.0},{3,4,-1.0},{4,3,-1.0},{4,4,4.0}
    };
    const std::vector<SparseTriplet> m_entries{
        {0,0,2.0},{1,1,2.0},{2,2,2.0},{3,3,2.0},{4,4,2.0}
    };
    const auto k = SparseMatrix::from_triplets(5, k_entries);
    const auto m = SparseMatrix::from_triplets(5, m_entries);
    const auto reference = solve_lowest_sparse_modes(k, m, 2, 500, 1e-10, false);
    require(reference.converged, "unpermuted reference solve did not converge");
    require(reference.vectors.size() >= 2, "reference returned too few vectors");

    const auto seeded = solve_lowest_sparse_modes(
        k, m, 2, 500, 1e-10, true, reference.vectors);
    require(seeded.converged, "RCM seeded solve did not converge");
    require(seeded.vectors.size() >= 2, "RCM seeded solve returned too few vectors");
    require(std::isfinite(seeded.final_relative_residual), "RCM seeded residual is not finite");
    require(seeded.final_relative_residual < 1e-7, "RCM seeded residual is too large");

    std::cout << "rcm_initial_vector_contract: PASS\n";
    std::cout << "reference_lambda0=" << reference.eigenvalues[0]
              << " seeded_lambda0=" << seeded.eigenvalues[0] << "\n";
    return 0;
}
