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

    // Deterministic 30-case invariance sweep. Each candidate spans exactly
    // the same plane, but its basis is independently rotated, permuted,
    // sign-flipped and rescaled. The metric subspace score must not depend on
    // the eigensolver's arbitrary basis choice inside a degenerate pair.
    constexpr int invariance_cases = 30;
    double worst_invariance_error = 0.0;
    for (int k = 0; k < invariance_cases; ++k) {
        const double angle = 0.173 * static_cast<double>(k + 1);
        const double cs = std::cos(angle);
        const double sn = std::sin(angle);
        const double scale0 = (k % 2 == 0 ? 1.0 : -1.0) * (0.2 + 0.13 * k);
        const double scale1 = (k % 3 == 0 ? -1.0 : 1.0) * (0.4 + 0.07 * k);
        std::vector<std::vector<double>> transformed{
            {scale0 * cs, scale0 * sn, 0.0},
            {-scale1 * sn, scale1 * cs, 0.0}};
        if (k % 2 != 0) std::swap(transformed[0], transformed[1]);
        const double score = mass_weighted_subspace_similarity(reference, transformed, mass);
        const double case_error = std::abs(score - 1.0);
        worst_invariance_error = std::max(worst_invariance_error, case_error);
        if (!(case_error < 1.0e-11) || !std::isfinite(score)) {
            std::cerr << std::setprecision(17)
                      << "basis invariance case " << k << " failed: score=" << score
                      << " error=" << case_error << '\n';
            return 1;
        }
    }

    // A controlled physical subspace change has a known minimum principal
    // overlap cos(theta). Sweep 30 angles to verify ranking and calibration.
    double worst_angle_error = 0.0;
    double previous_score = 1.0;
    for (int k = 1; k <= invariance_cases; ++k) {
        const double theta = 0.0025 * static_cast<double>(k);
        const double ct = std::cos(theta);
        const double st = std::sin(theta);
        const std::vector<std::vector<double>> tilted{
            {ct, 0.0, st},
            {0.0, 1.0, 0.0}};
        const double score = mass_weighted_subspace_similarity(reference, tilted, mass);
        const double expected_score = ct;
        worst_angle_error = std::max(worst_angle_error, std::abs(score - expected_score));
        if (!(std::abs(score - expected_score) < 1.0e-10) || score > previous_score + 1.0e-12) {
            std::cerr << std::setprecision(17)
                      << "principal-angle sweep case " << k << " failed: expected="
                      << expected_score << " measured=" << score << '\n';
            return 1;
        }
        previous_score = score;
    }

    std::cout << "modal_subspace_similarity_contract: PASS\n"
              << std::setprecision(17) << "expected_min_principal_overlap=" << expected
              << " measured=" << similarity
              << " error=" << error
              << " self=" << self
              << " invariance_cases=" << invariance_cases
              << " worst_invariance_error=" << worst_invariance_error
              << " angle_cases=" << invariance_cases
              << " worst_angle_error=" << worst_angle_error << '\n';
    return 0;
}
