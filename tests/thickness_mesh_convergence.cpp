#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
struct Level {
    int thickness_layers = 0;
    size_t nodes = 0;
    size_t elements = 0;
    BellMeshData mesh;
    BellOptimizer::ModalData modes;
};

bool makeLevel(const BellMeshData& base, int profile_nodes,
               int layers, const BellMaterial& material, Level& out) {
    if (layers < 1 || base.outer_profile.size() != static_cast<size_t>(profile_nodes) ||
        base.inner_profile.size() != static_cast<size_t>(profile_nodes)) return false;
    out.mesh.outer_profile = base.outer_profile;
    out.mesh.inner_profile = base.inner_profile;
    out.mesh.total_arc_length = base.total_arc_length;
    out.mesh.min_thickness = base.min_thickness;
    out.mesh.max_thickness = base.max_thickness;
    out.mesh.nodes.reserve(static_cast<size_t>(profile_nodes) * static_cast<size_t>(layers + 1));
    out.mesh.elements.reserve(static_cast<size_t>(profile_nodes - 1) * static_cast<size_t>(layers));

    // Subdivide the existing wall cross-section without changing either boundary.
    // Every intermediate node lies on the straight line between corresponding
    // outer/inner boundary nodes; only the FEM interpolation resolution changes.
    for (int layer = 0; layer <= layers; ++layer) {
        const double a = static_cast<double>(layer) / static_cast<double>(layers);
        for (int i = 0; i < profile_nodes; ++i) {
            const auto& o = base.outer_profile[static_cast<size_t>(i)];
            const auto& in = base.inner_profile[static_cast<size_t>(i)];
            Node node;
            node.id = static_cast<int>(out.mesh.nodes.size()) + 1;
            node.r = o.r + a * (in.r - o.r);
            node.z = o.z + a * (in.z - o.z);
            out.mesh.nodes.push_back(node);
        }
    }
    for (int layer = 0; layer < layers; ++layer) {
        const int row_outer = layer * profile_nodes;
        const int row_inner = (layer + 1) * profile_nodes;
        for (int i = 0; i < profile_nodes - 1; ++i) {
            Element e;
            e.n1 = row_inner + i + 1;
            e.n2 = row_outer + i + 1;
            e.n3 = row_outer + i + 2;
            e.n4 = row_inner + i + 2;
            out.mesh.elements.push_back(e);
        }
    }
    out.thickness_layers = layers;
    out.nodes = out.mesh.nodes.size();
    out.elements = out.mesh.elements.size();
    BellOptimizer::InternalModalSolver solver(true, 10, true);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{out.mesh, material, ".", "thickness_convergence", &modal_influence_1};
    return solver.solve(problem, out.modes) && out.modes.modes.size() >= 10;
}

struct ShapeMetrics { double correlation = 0.0; double error = 1.0; };
ShapeMetrics compareShape(const Level& a, const Level& b, int profile_nodes, size_t mode_index) {
    const auto& ma = a.modes.modes[mode_index];
    const auto& mb = b.modes.modes[mode_index];
    double dot = 0.0, na = 0.0, nb = 0.0;
    std::vector<double> va, vb;
    // Compare both fields at the same profile stations and at five normalized
    // through-thickness positions. Linear sampling is exact for nodal positions
    // and is used only as a diagnostic correspondence, not as an FE solution.
    for (int k = 0; k <= 4; ++k) {
        const double t = static_cast<double>(k) / 4.0;
        for (int i = 0; i < profile_nodes; ++i) {
            auto sample = [&](const Level& level, const BellOptimizer::ModalMode& mode,
                              const std::vector<double>& field) {
                const double x = t * level.thickness_layers;
                const int j0 = std::min(level.thickness_layers - 1, static_cast<int>(std::floor(x)));
                const int j1 = std::min(level.thickness_layers, j0 + 1);
                const double alpha = t >= 1.0 ? 1.0 : x - j0;
                const size_t id0 = static_cast<size_t>(j0 * profile_nodes + i);
                const size_t id1 = static_cast<size_t>(j1 * profile_nodes + i);
                (void)mode;
                return field[id0] + alpha * (field[id1] - field[id0]);
            };
            va.push_back(sample(a, ma, ma.nodal_radial));
            va.push_back(sample(a, ma, ma.nodal_axial));
            vb.push_back(sample(b, mb, mb.nodal_radial));
            vb.push_back(sample(b, mb, mb.nodal_axial));
        }
    }
    for (size_t i = 0; i < va.size(); ++i) {
        dot += va[i] * vb[i]; na += va[i] * va[i]; nb += vb[i] * vb[i];
    }
    if (!(na > 0.0) || !(nb > 0.0)) return {};
    const double scale = dot / na;
    double residual = 0.0;
    for (size_t i = 0; i < va.size(); ++i) {
        const double d = vb[i] - scale * va[i]; residual += d * d;
    }
    return {std::abs(dot) / std::sqrt(na * nb), std::sqrt(residual / nb)};
}
}

int main(int argc, char** argv) {
    const std::filesystem::path out_path = argc > 1 ? argv[1] : "thickness_mesh_convergence.csv";
    constexpr int profile_nodes = 31;
    const std::vector<int> layer_counts{1, 2, 3, 4};
    BellMaterial material;
    BellMeshConfig config; config.nNodes = profile_nodes;
    BellGenerator generator(material, config);
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) {
        std::cerr << "Failed to build fixed corrected geometry\n"; return 1;
    }
    const auto base = generator.getMeshData();
    std::vector<Level> levels;
    for (int layers : layer_counts) {
        Level level;
        if (!makeLevel(base, profile_nodes, layers, material, level)) {
            std::cerr << "Thickness-layer solve failed: layers=" << layers << '\n'; return 2;
        }
        std::cout << "profile_nodes=" << profile_nodes << ", thickness_layers=" << layers
                  << ", nodes=" << level.nodes << ", elements=" << level.elements
                  << ", f1=" << std::setprecision(10)
                  << level.modes.modes.front().frequency_hz << " Hz\n";
        levels.push_back(std::move(level));
    }

    std::ofstream csv(out_path);
    if (!csv) { std::cerr << "Cannot write " << out_path << '\n'; return 3; }
    csv << "profile_nodes,coarse_thickness_layers,fine_thickness_layers,mode,coarse_hz,fine_hz,relative_frequency_change,absolute_shape_correlation,normalized_shape_l2_error\n";
    csv << std::setprecision(12);
    double global_max_delta = 0.0, global_max_shape_error = 0.0, global_min_corr = 1.0;
    for (size_t l = 1; l < levels.size(); ++l) {
        const auto& coarse = levels[l - 1]; const auto& fine = levels[l];
        for (size_t m = 0; m < 10; ++m) {
            const double fc = coarse.modes.modes[m].frequency_hz;
            const double ff = fine.modes.modes[m].frequency_hz;
            const double delta = std::abs(ff - fc) / std::abs(ff);
            const auto shape = compareShape(coarse, fine, profile_nodes, m);
            global_max_delta = std::max(global_max_delta, delta);
            global_max_shape_error = std::max(global_max_shape_error, shape.error);
            global_min_corr = std::min(global_min_corr, shape.correlation);
            csv << profile_nodes << ',' << coarse.thickness_layers << ',' << fine.thickness_layers << ','
                << m + 1 << ',' << fc << ',' << ff << ',' << delta << ','
                << shape.correlation << ',' << shape.error << '\n';
        }
    }
    const auto summary_path = out_path.parent_path() / "thickness_mesh_convergence_summary.txt";
    std::ofstream summary(summary_path);
    if (!summary) { std::cerr << "Cannot write " << summary_path << '\n'; return 4; }
    summary << std::setprecision(12)
            << "geometry=corrected_inner_normal_offset\n"
            << "profile_nodes_fixed=" << profile_nodes << '\n'
            << "thickness_layers_sequence=1,2,3,4\n"
            << "boundary_geometry=unchanged; internal_nodes_linearly_subdivide_corresponding_wall_segments\n"
            << "modes_compared=10\n"
            << "fundamental_frequency_hz_by_thickness_layers\n";
    for (const auto& level : levels)
        summary << level.thickness_layers << ',' << level.modes.modes.front().frequency_hz << '\n';
    const auto& coarse = levels[levels.size() - 2]; const auto& fine = levels.back();
    double last_max_delta = 0.0, last_max_shape_error = 0.0, last_min_corr = 1.0;
    for (size_t m = 0; m < 10; ++m) {
        const double fc = coarse.modes.modes[m].frequency_hz;
        const double ff = fine.modes.modes[m].frequency_hz;
        last_max_delta = std::max(last_max_delta, std::abs(ff - fc) / std::abs(ff));
        const auto shape = compareShape(coarse, fine, profile_nodes, m);
        last_max_shape_error = std::max(last_max_shape_error, shape.error);
        last_min_corr = std::min(last_min_corr, shape.correlation);
    }
    summary << "\nfinal_adjacent_refinement=3_to_4_layers\n"
            << "max_frequency_change_modes_1_to_10_percent=" << last_max_delta * 100.0 << '\n'
            << "min_absolute_shape_correlation_modes_1_to_10=" << last_min_corr << '\n'
            << "max_normalized_shape_l2_error_modes_1_to_10_percent=" << last_max_shape_error * 100.0 << '\n'
            << "all_adjacent_pairs_max_frequency_change_percent=" << global_max_delta * 100.0 << '\n'
            << "all_adjacent_pairs_min_absolute_shape_correlation=" << global_min_corr << '\n'
            << "all_adjacent_pairs_max_normalized_shape_l2_error_percent=" << global_max_shape_error * 100.0 << '\n'
            << "limitation=profile_resolution_fixed_at_31; this isolates thickness refinement but does_not_establish_full_mesh_independence\n";
    std::cout << "Max adjacent frequency change, modes 1-10: " << global_max_delta * 100.0 << "%\n"
              << "Min adjacent absolute shape correlation: " << global_min_corr << '\n'
              << "Max adjacent normalized shape L2 error: " << global_max_shape_error * 100.0 << "%\n"
              << "Results: " << out_path << '\n';
    return csv && summary ? 0 : 5;
}
