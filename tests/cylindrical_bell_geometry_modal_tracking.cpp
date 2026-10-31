#include "bell_generator.h"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/cylindrical_modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

// Build a numerically useful Q4 meridional mesh from the actual BellGeometry
// profile. BellGenerator's production mesh is intentionally coarse in the
// thickness direction; this audit mesh samples the same inner/outer profiles
// with several radial layers so the vector FEM has enough resolution to carry
// non-axisymmetric families without changing the underlying bell geometry.
static BellMeshData make_refined_profile_mesh(const BellGenerator& generator,
                                              int radial_layers, int axial_layers) {
    BellMeshData mesh;
    for (int j = 0; j <= axial_layers; ++j) {
        const double z0 = 0.05 * generator.getGeometry().H;
        const double z1 = generator.getGeometry().H;
        const double z = z0 + (z1 - z0) * static_cast<double>(j) /
                         static_cast<double>(axial_layers);
        const double inner = generator.getInnerRadius(z);
        const double outer = generator.getOuterRadius(z);
        require(std::isfinite(inner) && std::isfinite(outer) && outer >= inner,
                "invalid real-bell profile sample");
        for (int i = 0; i <= radial_layers; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(radial_layers);
            const double r = inner + t * (outer - inner);
            mesh.nodes.push_back({static_cast<int>(mesh.nodes.size()) + 1, r * 1e-3, z * 1e-3});
        }
    }
    const int stride = radial_layers + 1;
    for (int j = 0; j < axial_layers; ++j) {
        for (int i = 0; i < radial_layers; ++i) {
            const int n1 = j * stride + i + 1;
            const int n2 = n1 + 1;
            const int n4 = (j + 1) * stride + i + 1;
            const int n3 = n4 + 1;
            mesh.elements.push_back({n1, n2, n3, n4});
        }
    }
    return mesh;
}

static std::vector<int> clamp_lower_edge(const BellMeshData& mesh, std::size_t dpn) {
    double z_min = mesh.nodes.front().z;
    for (const auto& node : mesh.nodes) z_min = std::min(z_min, node.z);
    std::vector<int> fixed;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (std::abs(mesh.nodes[i].z - z_min) < 1e-12) {
            for (std::size_t d = 0; d < dpn; ++d)
                fixed.push_back(static_cast<int>(i * dpn + d));
        }
    }
    return fixed;
}

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;

    // This is a real production BellGeometry used by the circular-limit audit.
    // No synthetic radii are introduced; only the meridional mesh resolution is
    // increased for this numerical family-tracking contract.
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    require(generator.validateGeometry().passed && generator.buildMesh(),
            "real bell geometry generation failed");

    const BellMeshData mesh = make_refined_profile_mesh(generator, 3, 12);
    require(mesh.nodes.size() == 52 && mesh.elements.size() == 36,
            "unexpected refined real-bell mesh dimensions");

    const double E_pa = 1.0e11;
    const double rho_kg_m3 = 7800.0;
    const std::vector<int> harmonics{0, 1};
    constexpr int requested = 4;

    const auto reference_system = assemble_cylindrical_global_system(
        mesh, E_pa, material.nu, rho_kg_m3, harmonics, 64, 0.0, 1.0);
    const auto fixed = clamp_lower_edge(mesh, reference_system.dofs_per_node);
    const auto reference = solve_cylindrical_global_modes(
        reference_system, fixed, requested, 5000, 1e-2, true);
    std::cout << "reference_count=" << reference.frequencies_hz.size()
              << " converged=" << reference.converged
              << " residual=" << reference.relative_residual << "\n";
    require(reference.converged, "real-bell q=1 vector solve did not converge");
    require(reference.frequencies_hz.size() >= static_cast<std::size_t>(requested),
            "real-bell q=1 returned too few modes");
    require(reference.relative_residual < 2e-4, "real-bell q=1 residual too large");

    std::cout << "real_bell_q1_modes:";
    for (double f : reference.frequencies_hz) std::cout << ' ' << f;
    std::cout << "\n";

    // Find a genuine circular-limit degeneracy rather than assuming its rank.
    std::size_t pair_a = 0, pair_b = 1;
    double pair_split = 1e9;
    for (std::size_t i = 0; i + 1 < reference.frequencies_hz.size(); ++i) {
        const double split = std::abs(reference.frequencies_hz[i + 1] - reference.frequencies_hz[i]) /
                             std::max(reference.frequencies_hz[i], 1e-12);
        if (split < pair_split) {
            pair_split = split;
            pair_a = i;
            pair_b = i + 1;
        }
    }
    require(pair_split < 2e-3, "real-bell q=1 has no sufficiently degenerate modal pair");

    const auto& metric_mass = reference_system.mass;
    require(mass_weighted_subspace_similarity(
                reference, {pair_a, pair_b}, reference, {pair_a, pair_b}, metric_mass) > 1.0 - 1e-9,
            "real-bell reference pair is not self-consistent");

    std::cout << "reference_pair=" << pair_a << ',' << pair_b
              << " relative_split=" << pair_split << "\n";

    double previous_similarity = 1.0;
    std::vector<std::vector<double>> previous_modes = reference.mode_shapes;
    for (double q : {0.95, 0.90}) {
        const auto system = assemble_cylindrical_global_system(
            mesh, E_pa, material.nu, rho_kg_m3, harmonics, 64, 0.0, q);
        const auto candidate = solve_cylindrical_global_modes(
            system, fixed, requested, 5000, 1e-2, true, previous_modes);
        require(candidate.converged, "real-bell elliptical vector solve did not converge");
        require(candidate.frequencies_hz.size() >= static_cast<std::size_t>(requested),
                "real-bell elliptical solve returned too few modes");
        require(candidate.relative_residual < 2e-4, "real-bell elliptical residual too large");

        double best_similarity = -1.0;
        std::size_t best_i = 0, best_j = 1;
        for (std::size_t i = 0; i < candidate.frequencies_hz.size(); ++i) {
            for (std::size_t j = i + 1; j < candidate.frequencies_hz.size(); ++j) {
                const double score = mass_weighted_subspace_similarity(
                    reference, {pair_a, pair_b}, candidate, {i, j}, metric_mass);
                if (score > best_similarity) {
                    best_similarity = score;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        const double ov00 = mass_weighted_modal_overlap(reference, pair_a, candidate, best_i, metric_mass);
        const double ov01 = mass_weighted_modal_overlap(reference, pair_a, candidate, best_j, metric_mass);
        const double ov10 = mass_weighted_modal_overlap(reference, pair_b, candidate, best_i, metric_mass);
        const double ov11 = mass_weighted_modal_overlap(reference, pair_b, candidate, best_j, metric_mass);
        const double individual = std::max(ov00 + ov11, ov01 + ov10) / 2.0;
        const double split = std::abs(candidate.frequencies_hz[best_j] - candidate.frequencies_hz[best_i]) /
                             std::max(std::min(candidate.frequencies_hz[best_i], candidate.frequencies_hz[best_j]), 1e-12);

        require(best_similarity > 0.45, "real-bell modal family lost continuity");
        require(best_similarity <= 1.0 + 1e-8, "real-bell subspace similarity exceeded unity");
        require(split > pair_split + 1e-4, "ellipticity did not split the tracked real-bell family");
        require(best_similarity > previous_similarity - 0.35,
                "real-bell modal-family overlap collapsed abruptly");
        previous_similarity = best_similarity;
        previous_modes = candidate.mode_shapes;

        std::cout << "q=" << q
                  << " tracked_pair=" << best_i << ',' << best_j
                  << " pair_freq_Hz=" << candidate.frequencies_hz[best_i]
                  << ',' << candidate.frequencies_hz[best_j]
                  << " split=" << split
                  << " subspace_similarity=" << best_similarity
                  << " individual_overlap=" << individual
                  << " residual=" << candidate.relative_residual << "\n";
    }

    std::cout << "cylindrical_bell_geometry_modal_tracking: PASS\n";
    return 0;
}
