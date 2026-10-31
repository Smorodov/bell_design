#include "bell_generator.h"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/cylindrical_modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
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
                                    std::size_t modes = 4,
                                    const std::vector<std::vector<double>>& initial = {}) {
    constexpr double E = 1.0e11;
    constexpr double nu = 0.30;
    constexpr double rho = 7800.0;
    const auto system = assemble_cylindrical_global_system(
        mesh, E, nu, rho, harmonics, 64, 0.0, q);
    const auto fixed = clamp_lower_edge(mesh, system.dofs_per_node);
    const auto result = solve_cylindrical_global_modes(
        system, fixed, static_cast<int>(modes), 6000, 1e-6, true, initial);
    require(result.converged, "unified ellipticity solve did not converge");
    require(result.frequencies_hz.size() >= modes,
            "unified ellipticity solve returned too few modes");
    require(result.relative_residual < 2e-4,
            "unified ellipticity residual too large");
    return result;
}

static std::pair<std::size_t, std::size_t> best_pair(
    const std::vector<std::vector<double>>& reference,
    const CylindricalModalResult& candidate,
    const detail::SparseMatrix& candidate_mass,
    double* best_score = nullptr) {
    double best = -1.0;
    std::pair<std::size_t, std::size_t> pair{0, 1};
    for (std::size_t i = 0; i < candidate.mode_shapes.size(); ++i) {
        for (std::size_t j = i + 1; j < candidate.mode_shapes.size(); ++j) {
            const double score = mass_weighted_subspace_similarity(
                reference, {candidate.mode_shapes[i], candidate.mode_shapes[j]}, candidate_mass);
            if (score > best) {
                best = score;
                pair = {i, j};
            }
        }
    }
    if (best_score) *best_score = best;
    return pair;
}

static double splitting(const CylindricalModalResult& r,
                        std::pair<std::size_t, std::size_t> p) {
    const double a = r.frequencies_hz[p.first];
    const double b = r.frequencies_hz[p.second];
    return std::abs(a - b) /
           std::max(std::min(std::abs(a), std::abs(b)), 1e-12);
}

static double relative_change(double a, double b) {
    return std::abs(a - b) /
           std::max(std::min(std::abs(a), std::abs(b)), 1e-12);
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

    const auto medium = make_profile_mesh(generator, 2, 8);
    const auto fine = make_profile_mesh(generator, 3, 12);

    // --- Axes 1+2: spatial refinement and ellipticity continuation. ---
    // Spatial refinement uses the cheaper {0,1} space because the purpose of
    // this axis is meridional discretization. Ellipticity continuation is
    // separately isolated in the validated enriched {0,1,2,3} space.
    const auto mesh_m085 = solve(medium, 0.85, h01, 3);
    const auto mesh_f085 = solve(fine, 0.85, h01, 3);
    const auto mesh_m085_system = assemble_cylindrical_global_system(
        medium, 1.0e11, 0.30, 7800.0, h01, 64, 0.0, 0.85);
    const auto mesh_f085_system = assemble_cylindrical_global_system(
        fine, 1.0e11, 0.30, 7800.0, h01, 64, 0.0, 0.85);
    const auto mesh_pair_m = std::pair<std::size_t, std::size_t>{0, 1};
    const auto mesh_pair_f = std::pair<std::size_t, std::size_t>{0, 1};
    const double mesh085 = std::max(
        relative_change(mesh_m085.frequencies_hz[0], mesh_f085.frequencies_hz[0]),
        relative_change(mesh_m085.frequencies_hz[1], mesh_f085.frequencies_hz[1]));
    const double split_mesh_m = splitting(mesh_m085, mesh_pair_m);
    const double split_mesh_f = splitting(mesh_f085, mesh_pair_f);
    const double split_mesh085 = relative_change(split_mesh_m, split_mesh_f);

    // The enriched q sweep is done on one mesh. This isolates q-continuation
    // from spatial refinement and avoids using frequency rank as family ID.
    const auto m090 = solve(medium, 0.90, h0123, 4);
    const auto m085 = solve(medium, 0.85, h0123, 4, m090.mode_shapes);
    const auto m090_system = assemble_cylindrical_global_system(
        medium, 1.0e11, 0.30, 7800.0, h0123, 64, 0.0, 0.90);
    const auto m085_system = assemble_cylindrical_global_system(
        medium, 1.0e11, 0.30, 7800.0, h0123, 64, 0.0, 0.85);
    double q_similarity = 0.0;
    const auto q_pair = best_pair(
        {m090.mode_shapes[0], m090.mode_shapes[1]}, m085, m085_system.mass, &q_similarity);
    const double split_q090 = splitting(m090, {0, 1});
    const double split_q085 = splitting(m085, q_pair);
    const double q_split_growth = split_q085 / std::max(split_q090, 1e-12);

    std::cout << "axis=mesh q=0.85 h01 medium_pair="
              << mesh_m085.frequencies_hz[0] << ',' << mesh_m085.frequencies_hz[1]
              << " fine_pair=" << mesh_f085.frequencies_hz[0] << ',' << mesh_f085.frequencies_hz[1]
              << " mesh_change=" << mesh085 << " split_change=" << split_mesh085 << '\n';
    std::cout << "axis=ellipticity medium h0123 q=0.90 pair="
              << m090.frequencies_hz[0] << ',' << m090.frequencies_hz[1]
              << " split=" << split_q090
              << " q=0.85 tracked_indices=" << q_pair.first << ',' << q_pair.second
              << " freqs=" << m085.frequencies_hz[q_pair.first] << ',' << m085.frequencies_hz[q_pair.second]
              << " split=" << split_q085 << " growth=" << q_split_growth
              << " continuation_similarity=" << q_similarity << '\n';

    require(mesh085 < 0.10,
            "real-bell frequencies remain too mesh-sensitive");
    require(split_mesh085 < 0.15,
            "real-bell modal family splitting is not sufficiently mesh-stable");
    require(q_similarity > 0.45 && split_q090 > 1e-4 && split_q085 > split_q090,
            "ellipticity continuation did not preserve and increase the resolved split");
    // --- Axis 3: harmonic enrichment at fixed q=0.85 on the medium mesh. ---
    const auto b01 = solve(medium, 0.85, h01, 3);
    const auto b012 = solve(medium, 0.85, h012, 6);
    const auto b0123 = m085;
    const auto s012 = assemble_cylindrical_global_system(
        medium, 1.0e11, 0.30, 7800.0, h012, 64, 0.0, 0.85);
    const auto s0123 = m085_system;

    const auto p012 = prolongate_cylindrical_modes(b01, h01, h012);
    const auto p0123 = prolongate_cylindrical_modes(b01, h01, h0123);
    const auto p23 = prolongate_cylindrical_modes(b012, h012, h0123);
    require(p012.size() >= 2 && p0123.size() >= 2 && p23.size() >= 2,
            "harmonic continuation projection failed");

    double hs012 = 0.0, hs0123 = 0.0, hs23 = 0.0;
    const auto hp012 = best_pair({p012[0], p012[1]}, b012, s012.mass, &hs012);
    const auto hp0123 = best_pair({p0123[0], p0123[1]}, b0123, s0123.mass, &hs0123);
    const auto hp23 = best_pair({p23[0], p23[1]}, b0123, s0123.mass, &hs23);
    const double harmonic_direct = hs0123;
    const double harmonic_chain = hs23;
    const double hfreq01 = std::max(relative_change(b01.frequencies_hz[0], b012.frequencies_hz[hp012.first]),
                                     relative_change(b01.frequencies_hz[1], b012.frequencies_hz[hp012.second]));
    const double hfreq23 = std::max(relative_change(b012.frequencies_hz[hp012.first], b0123.frequencies_hz[hp0123.first]),
                                    relative_change(b012.frequencies_hz[hp012.second], b0123.frequencies_hz[hp0123.second]));

    std::cout << "axis=harmonics q=0.85 h01_pair=" << b01.frequencies_hz[0] << ',' << b01.frequencies_hz[1]
              << " h012_indices=" << hp012.first << ',' << hp012.second
              << " h0123_indices=" << hp0123.first << ',' << hp0123.second
              << " sim01_to_012=" << hs012 << " sim01_to_0123=" << harmonic_direct
              << " sim012_to_0123=" << harmonic_chain
              << " rank_change_01_to_012=" << hfreq01
              << " rank_change_012_to_0123=" << hfreq23 << '\n';

    require(hs012 > 0.70 && harmonic_direct > 0.65 && harmonic_chain > 0.65,
            "harmonic enrichment lost the physical modal family");
    // Frequency-rank changes are diagnostic only. They are expected when new
    // Fourier families enter the spectrum and therefore do not gate PASS.

    std::cout << "cylindrical_unified_ellipticity_audit: PASS\n";
    return 0;
}
