#include "fem/cylindrical_modal_tracking.hpp"
#include <iomanip>

#include <cmath>
#include <iostream>
#include <vector>

using namespace BellOptimizer::FEM;
using BellOptimizer::detail::SparseMatrix;

int main() {
    // Identity mass metric. The reference subspace is span(e1,e2); the
    // candidate is span(e1, cos(theta)e2 + sin(theta)e3). Their principal
    // overlaps are exactly 1 and cos(theta). Choose a very small separation so
    // the smaller singular value is close to one, where the naive
    // trace-minus-discriminant formula suffers catastrophic cancellation.
    const double eps = 1.0e-8;
    const double c = 1.0 - eps;
    const double s = std::sqrt(std::max(0.0, 1.0 - c * c));

    SparseMatrix mass = SparseMatrix::from_triplets(
        3, {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}});

    const std::vector<std::vector<double>> reference{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}};
    const std::vector<std::vector<double>> candidate{
        {1.0, 0.0, 0.0},
        {0.0, c, s}};

    const double similarity = mass_weighted_subspace_similarity(reference, candidate, mass);
    const double expected = c;
    const double error = std::abs(similarity - expected);

    // The stable implementation should retain the small 1e-8 separation to
    // substantially better than 1e-9 absolute error.
    if (!(error < 1.0e-10) || !(similarity < 1.0) || !std::isfinite(similarity)) {
        std::cerr << std::setprecision(17) << "modal_subspace_similarity_contract: FAIL\n"
                  << "expected=" << expected
                  << " similarity=" << similarity
                  << " error=" << error << '\n';
        return 1;
    }

    // Identical subspaces must remain exactly/near-exactly self-consistent.
    const double self = mass_weighted_subspace_similarity(reference, reference, mass);
    if (!(self > 1.0 - 1.0e-12) || !(self <= 1.0 + 1.0e-12)) {
        std::cerr << "modal_subspace_similarity_contract: FAIL self=" << self << '\n';
        return 1;
    }

    std::cout << "modal_subspace_similarity_contract: PASS\n"
              << std::setprecision(17) << "expected_min_principal_overlap=" << expected
              << " measured=" << similarity
              << " error=" << error
              << " self=" << self << '\n';
    return 0;
}
