#include "core/bell_types.h"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/cylindrical_modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

CylindricalModalResult solve(const CylindricalGlobalSystem& system,
                             const CylindricalConstraintTransform& transform,
                             bool rcm) {
    auto result = solve_cylindrical_global_modes(system, transform, 24, 2400, 1e-6, rcm);
    require(result.converged && result.frequencies_hz.size() >= 12,
            "full cylindrical FEM solve did not converge or returned too few modes");
    require(result.relative_residual < 1e-5, "full cylindrical FEM residual is too large");
    return result;
}

std::vector<std::pair<std::size_t, std::size_t>> find_degenerate_pairs(
    const CylindricalModalResult& result, std::size_t count, double relative_gap) {
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    const std::size_t n = std::min(count, result.frequencies_hz.size());
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const double a = result.frequencies_hz[i];
        const double b = result.frequencies_hz[i + 1];
        if (std::abs(b - a) / std::max({std::abs(a), std::abs(b), 1e-12}) <= relative_gap)
            pairs.emplace_back(i, i + 1);
    }
    return pairs;
}

std::vector<std::vector<double>> pair_vectors(const CylindricalModalResult& result,
                                               std::size_t i, std::size_t j) {
    return {result.mode_shapes.at(i), result.mode_shapes.at(j)};
}
} // namespace

int main() {
    try {
        BellGeometry geometry;
        geometry.H = 540.0;
        geometry.R_mouth = 300.0;
        geometry.t_top = 14.0;
        geometry.t_bottom = 34.0;
        geometry.r_knots = {0.0, 0.085, 0.19, 0.43, 1.0};
        geometry.ellipticity = 1.0;
        BellMaterial material;
        BellMeshConfig mesh;
        mesh.nNodes = 15;

        const auto model = make_cylindrical_bell_model(geometry, material, mesh);
        const std::vector<int> harmonics{0, 1, 2, 3};
        const auto system = assemble_cylindrical_bell_system(model, harmonics, 32);
        const auto transform = make_cylindrical_axis_regularity_transform(
            model.mesh_si, system.dofs_per_node, harmonics);

        // The two runs solve exactly the same physical generalized eigenproblem,
        // but use different sparse orderings. Circular-bell cosine/sine pairs
        // may return a different orthonormal basis inside the same eigenspace.
        const auto ordered = solve(system, transform, true);
        const auto natural = solve(system, transform, false);
        const auto pairs = find_degenerate_pairs(ordered, 24, 1e-6);
        require(!pairs.empty(), "no near-degenerate circular-bell modal pair found");

        std::cout << std::setprecision(12)
                  << "consistent_mass_dimension=" << system.mass.n
                  << " consistent_mass_nnz=" << system.mass.values.size()
                  << " reference_modes=" << ordered.frequencies_hz.size()
                  << " alternate_modes=" << natural.frequencies_hz.size()
                  << " detected_pairs=" << pairs.size() << '\n';

        double minimum_subspace = 1.0;
        double maximum_frequency_pair_error = 0.0;
        for (const auto& pair : pairs) {
            double best = -1.0;
            std::size_t best_i = 0, best_j = 1;
            for (std::size_t i = 0; i + 1 < natural.mode_shapes.size(); ++i) {
                for (std::size_t j = i + 1; j < natural.mode_shapes.size(); ++j) {
                    // Compare explicit two-mode families; candidate indices need not be adjacent.
                    const double candidate_similarity = mass_weighted_subspace_similarity(
                        pair_vectors(ordered, pair.first, pair.second),
                        pair_vectors(natural, i, j), system.mass);
                    if (candidate_similarity > best) {
                        best = candidate_similarity;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            const double individual00 = mass_weighted_modal_overlap(
                ordered, pair.first, natural, best_i, system.mass);
            const double individual01 = mass_weighted_modal_overlap(
                ordered, pair.first, natural, best_j, system.mass);
            const double individual10 = mass_weighted_modal_overlap(
                ordered, pair.second, natural, best_i, system.mass);
            const double individual11 = mass_weighted_modal_overlap(
                ordered, pair.second, natural, best_j, system.mass);
            const double individual_best = std::max(
                std::min(individual00, individual11), std::min(individual01, individual10));
            const double frequency_error = std::max(
                std::abs(ordered.frequencies_hz[pair.first] - natural.frequencies_hz[best_i]) /
                    std::max(ordered.frequencies_hz[pair.first], 1e-12),
                std::abs(ordered.frequencies_hz[pair.second] - natural.frequencies_hz[best_j]) /
                    std::max(ordered.frequencies_hz[pair.second], 1e-12));
            minimum_subspace = std::min(minimum_subspace, best);
            maximum_frequency_pair_error = std::max(maximum_frequency_pair_error, frequency_error);
            std::cout << "pair=" << pair.first << ',' << pair.second
                      << " frequencies_hz=" << ordered.frequencies_hz[pair.first] << ','
                      << ordered.frequencies_hz[pair.second]
                      << " matched_pair=" << best_i << ',' << best_j
                      << " min_principal_overlap=" << best
                      << " best_individual_assignment_overlap=" << individual_best
                      << " max_pair_frequency_relative_error=" << frequency_error << '\n';
        }
        std::cout << "minimum_rank2_subspace_similarity=" << minimum_subspace
                  << " maximum_pair_frequency_relative_error=" << maximum_frequency_pair_error << '\n';
        require(minimum_subspace > 0.999,
                "consistent-mass rank-2 modal subspace did not remain invariant across solver orderings");
        std::cout << "round_bell_consistent_mass_subspace_probe: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "round_bell_consistent_mass_subspace_probe: FAIL: " << error.what() << '\n';
        return 1;
    }
}
