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

static BellMeshData make_profile_mesh(const BellGenerator& generator,
                                      int radial_layers, int axial_layers) {
    BellMeshData mesh;
    const auto& g = generator.getGeometry();
    const double z0 = 0.05 * g.H;
    for (int j = 0; j <= axial_layers; ++j) {
        const double z = z0 + (g.H - z0) * static_cast<double>(j) /
                         static_cast<double>(axial_layers);
        const double inner = generator.getInnerRadius(z);
        const double outer = generator.getOuterRadius(z);
        require(std::isfinite(inner) && std::isfinite(outer) && outer >= inner,
                "invalid real-bell profile sample");
        for (int i = 0; i <= radial_layers; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(radial_layers);
            mesh.nodes.push_back({static_cast<int>(mesh.nodes.size()) + 1,
                                  (inner + t * (outer - inner)) * 1e-3, z * 1e-3});
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
        if (std::abs(mesh.nodes[i].z - z_min) < 1e-12)
            for (std::size_t d = 0; d < dpn; ++d)
                fixed.push_back(static_cast<int>(i * dpn + d));
    }
    return fixed;
}

static BellGeometry reference_geometry() {
    BellGeometry g;
    g.H = 658.3175709034034;
    g.R_mouth = 355.9758748202059;
    g.t_top = 22.85675347179014;
    g.t_bottom = 46.62524667849419;
    g.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    return g;
}

static CylindricalModalResult solve(const BellMeshData& mesh, double q,
                                    const std::vector<int>& harmonics,
                                    std::size_t modes,
                                    const std::vector<std::vector<double>>& initial = {}) {
    constexpr double E = 1.0e11;
    constexpr double nu = 0.30;
    constexpr double rho = 7800.0;
    const auto system = assemble_cylindrical_global_system(
        mesh, E, nu, rho, harmonics, 64, 0.0, q);
    const auto fixed = clamp_lower_edge(mesh, system.dofs_per_node);
    const auto result = solve_cylindrical_global_modes(
        system, fixed, static_cast<int>(modes), 6000, 1e-6, true, initial);
    require(result.converged, "harmonic continuation solve did not converge");
    require(result.frequencies_hz.size() >= modes, "harmonic continuation returned too few modes");
    require(result.relative_residual < 2e-4, "harmonic continuation residual too large");
    return result;
}

static std::pair<std::size_t, std::size_t> best_pair(
    const std::vector<std::vector<double>>& projected,
    const CylindricalModalResult& candidate,
    const detail::SparseMatrix& candidate_mass) {
    double best = -1.0;
    std::pair<std::size_t, std::size_t> pair{0, 1};
    for (std::size_t i = 0; i < candidate.mode_shapes.size(); ++i) {
        for (std::size_t j = i + 1; j < candidate.mode_shapes.size(); ++j) {
            const double score = mass_weighted_subspace_similarity(
                projected, {candidate.mode_shapes[i], candidate.mode_shapes[j]}, candidate_mass);
            if (score > best) {
                best = score;
                pair = {i, j};
            }
        }
    }
    require(best > 0.35, "harmonic modal family lost continuity after enrichment");
    return pair;
}

static double splitting(const CylindricalModalResult& r, std::pair<std::size_t, std::size_t> p) {
    const double a = r.frequencies_hz[p.first];
    const double b = r.frequencies_hz[p.second];
    return std::abs(a - b) / std::max(std::min(std::abs(a), std::abs(b)), 1e-12);
}

int main() {
    const std::vector<int> h01{0, 1};
    const std::vector<int> h012{0, 1, 2};
    const std::vector<int> h0123{0, 1, 2, 3};

    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(reference_geometry());
    require(generator.validateGeometry().passed && generator.buildMesh(),
            "real bell geometry generation failed");
    const auto mesh = make_profile_mesh(generator, 2, 8);

    constexpr double q = 0.85;
    const auto base = solve(mesh, q, h01, 3);
    const auto enriched2 = solve(mesh, q, h012, 6);
    const auto enriched3 = solve(mesh, q, h0123, 8);

    // Exact nested-space embedding: {0,1} is a coordinate subspace of both
    // enriched bases. The candidate mass is the correct metric for the larger
    // physical space because ellipticity couples its Fourier harmonics.
    const auto projected2 = prolongate_cylindrical_modes(base, h01, h012);
    const auto projected3 = prolongate_cylindrical_modes(base, h01, h0123);
    const auto projected23 = prolongate_cylindrical_modes(enriched2, h012, h0123);
    require(projected2.size() >= 2 && projected3.size() >= 2 && projected23.size() >= 2,
            "failed to prolongate reference modal family");

    const auto system2 = assemble_cylindrical_global_system(
        mesh, 1.0e11, 0.30, 7800.0, h012, 64, 0.0, q);
    const auto system3 = assemble_cylindrical_global_system(
        mesh, 1.0e11, 0.30, 7800.0, h0123, 64, 0.0, q);
    const auto pair2 = best_pair({projected2[0], projected2[1]}, enriched2, system2.mass);
    const auto pair3 = best_pair({projected3[0], projected3[1]}, enriched3, system3.mass);
    const auto pair3_chain = best_pair({projected23[0], projected23[1]}, enriched3, system3.mass);
    const double sim2 = mass_weighted_subspace_similarity(
        {projected2[0], projected2[1]}, {enriched2.mode_shapes[pair2.first], enriched2.mode_shapes[pair2.second]}, system2.mass);
    const double sim3 = mass_weighted_subspace_similarity(
        {projected3[0], projected3[1]}, {enriched3.mode_shapes[pair3.first], enriched3.mode_shapes[pair3.second]}, system3.mass);
    const double chain23 = mass_weighted_subspace_similarity(
        {projected23[0], projected23[1]}, {enriched3.mode_shapes[pair3_chain.first], enriched3.mode_shapes[pair3_chain.second]}, system3.mass);
    const double split_base = splitting(base, {0, 1});
    const double split2 = splitting(enriched2, pair2);
    const double split3 = splitting(enriched3, pair3);

    std::cout << "q=" << q << '\n'
              << "h01_pair=" << base.frequencies_hz[0] << ',' << base.frequencies_hz[1]
              << " split=" << split_base << '\n'
              << "h012_pair_indices=" << pair2.first << ',' << pair2.second
              << " freqs=" << enriched2.frequencies_hz[pair2.first] << ',' << enriched2.frequencies_hz[pair2.second]
              << " similarity=" << sim2 << " split=" << split2 << '\n'
              << "h0123_pair_indices=" << pair3.first << ',' << pair3.second
              << " freqs=" << enriched3.frequencies_hz[pair3.first] << ',' << enriched3.frequencies_hz[pair3.second]
              << " similarity=" << sim3 << " chain_pair_indices=" << pair3_chain.first << ',' << pair3_chain.second
              << " chain_012_to_0123=" << chain23 << " split=" << split3 << '\n';

    // The point of this contract is family continuity, not equality of sorted
    // frequency ranks. A newly admitted harmonic can create lower modes and
    // reorder the spectrum. Require the same two-dimensional physical family
    // to remain recognizable in both enriched spaces.
    require(sim2 > 0.70, "{0,1}->{0,1,2} family similarity is too low");
    require(sim3 > 0.65, "{0,1}->{0,1,2,3} family similarity is too low");
    require(chain23 > 0.65, "{0,1,2}->{0,1,2,3} family similarity is too low");
    require(split2 > 1e-4 && split3 > 1e-4,
            "ellipticity did not preserve a resolvable tracked-family split");

    std::cout << "cylindrical_harmonic_family_continuation: PASS\n";
    return 0;
}
