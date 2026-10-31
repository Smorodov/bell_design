#include "fem/cylindrical_modal_tracking.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace BellOptimizer::FEM;
using BellOptimizer::detail::SparseMatrix;

int main() {
    const double delta = 1.0e-10;
    SparseMatrix mass = SparseMatrix::from_triplets(
        3, {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}});

    // The two input bases are individually almost linearly dependent but
    // span exactly the same plane.  The principal overlap must therefore
    // remain one.  This is a conditioning case: the second basis direction
    // is only delta away from the first and must not be discarded merely
    // because its residual is small in absolute terms.
    const std::vector<std::vector<double>> a{
        {1.0, 0.0, 0.0},
        {1.0, delta, 0.0}};
    const std::vector<std::vector<double>> b{
        {1.0, delta, 0.0},
        {1.0, 2.0 * delta, 0.0}};

    const double similarity = mass_weighted_subspace_similarity(a, b, mass);
    const double expected = 1.0;
    const double error = std::abs(similarity - expected);

    if (!std::isfinite(similarity) || error > 1.0e-8 || similarity < 0.0 || similarity > 1.0) {
        std::cerr << std::setprecision(17)
                  << "modal_subspace_conditioning_contract: FAIL\n"
                  << "expected=" << expected
                  << " measured=" << similarity
                  << " error=" << error << '\n';
        return 1;
    }

    const double self = mass_weighted_subspace_similarity(a, a, mass);
    if (!std::isfinite(self) || self < 1.0 - 1.0e-12 || self > 1.0 + 1.0e-12) {
        std::cerr << std::setprecision(17)
                  << "modal_subspace_conditioning_contract: FAIL self=" << self << '\n';
        return 1;
    }

    std::cout << "modal_subspace_conditioning_contract: PASS\n"
              << std::setprecision(17)
              << "expected_min_principal_overlap=" << expected
              << " measured=" << similarity
              << " error=" << error
              << " self=" << self << '\n';
    return 0;
}
