#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct Level {
    int profile_nodes = 0;
    size_t mesh_nodes = 0;
    size_t elements = 0;
    BellOptimizer::ModalData modes;
};

bool buildLevel(int n, const BellMaterial& material, Level& level) {
    BellMeshConfig config;
    config.nNodes = n;
    BellGenerator generator(material, config);
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    generator.setGeometry(geometry);
    const auto validation = generator.validateGeometry();
    if (!validation.passed) {
        std::cerr << "nNodes=" << n << ": geometry validation failed: " << validation.reason << '\n';
        return false;
    }
    if (!generator.buildMesh()) {
        std::cerr << "nNodes=" << n << ": mesh generation failed\n";
        return false;
    }
    const auto& mesh = generator.getMeshData();
    BellOptimizer::InternalModalSolver solver(true, 10, true);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{mesh, material, ".", "mesh_convergence", &modal_influence_1};
    if (!solver.solve(problem, level.modes) || level.modes.modes.size() < 10) {
        std::cerr << "nNodes=" << n << ": internal modal solve failed or returned fewer than 10 modes\n";
        return false;
    }
    level.profile_nodes = n;
    level.mesh_nodes = mesh.nodes.size();
    level.elements = mesh.elements.size();
    return true;
}

// Interpolate a mode's displacement samples by normalized profile index. The
// generator distributes profile samples by arc length, so this is a stable
// first-order correspondence for comparing different resolutions.
void sampleShape(const BellOptimizer::ModalMode& mode, int profile_nodes, double t,
                 double& ur_outer, double& uz_outer, double& ur_inner, double& uz_inner) {
    const double x = std::clamp(t, 0.0, 1.0) * static_cast<double>(profile_nodes - 1);
    const int i0 = std::min(profile_nodes - 1, static_cast<int>(std::floor(x)));
    const int i1 = std::min(profile_nodes - 1, i0 + 1);
    const double a = x - static_cast<double>(i0);
    auto interp = [&](const std::vector<double>& v, int offset) {
        const double p = v[static_cast<size_t>(offset + i0)];
        const double q = v[static_cast<size_t>(offset + i1)];
        return p + a * (q - p);
    };
    ur_outer = interp(mode.nodal_radial, 0);
    uz_outer = interp(mode.nodal_axial, 0);
    ur_inner = interp(mode.nodal_radial, profile_nodes);
    uz_inner = interp(mode.nodal_axial, profile_nodes);
}

struct ShapeMetrics { double correlation = 0.0; double error = 1.0; };
ShapeMetrics compareShape(const BellOptimizer::ModalMode& coarse, int nc,
                          const BellOptimizer::ModalMode& fine, int nf) {
    // Evaluate both fields at the fine profile's normalized sample positions.
    double dot = 0.0, norm_c = 0.0, norm_f = 0.0;
    std::vector<double> c, f;
    c.reserve(static_cast<size_t>(nf) * 4);
    f.reserve(static_cast<size_t>(nf) * 4);
    for (int i = 0; i < nf; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(nf - 1);
        double curo, cuzo, curi, cuzi, furo, fuzo, furi, fuzi;
        sampleShape(coarse, nc, t, curo, cuzo, curi, cuzi);
        sampleShape(fine, nf, t, furo, fuzo, furi, fuzi);
        c.insert(c.end(), {curo, cuzo, curi, cuzi});
        f.insert(f.end(), {furo, fuzo, furi, fuzi});
    }
    for (size_t i = 0; i < c.size(); ++i) {
        dot += c[i] * f[i]; norm_c += c[i] * c[i]; norm_f += f[i] * f[i];
    }
    if (!(norm_c > 0.0) || !(norm_f > 0.0)) return {};
    const double scale = dot / norm_c;
    double residual = 0.0;
    for (size_t i = 0; i < c.size(); ++i) {
        const double d = f[i] - scale * c[i]; residual += d * d;
    }
    return {std::abs(dot) / std::sqrt(norm_c * norm_f), std::sqrt(residual / norm_f)};
}
}

int main(int argc, char** argv) {
    const std::filesystem::path out_path = argc > 1 ? argv[1] : "mesh_convergence.csv";
    const std::vector<int> resolutions{41, 51, 61, 71, 81};
    BellMaterial material;
    std::vector<Level> levels;
    for (int n : resolutions) {
        Level level;
        if (!buildLevel(n, material, level)) return 1;
        std::cout << "nNodes=" << n << ", mesh nodes=" << level.mesh_nodes
                  << ", elements=" << level.elements << ", f1="
                  << std::setprecision(10) << level.modes.modes.front().frequency_hz << " Hz\n";
        levels.push_back(std::move(level));
    }

    std::ofstream csv(out_path);
    if (!csv) { std::cerr << "Cannot write " << out_path << '\n'; return 2; }
    csv << "coarse_profile_nodes,fine_profile_nodes,mode,coarse_hz,fine_hz,relative_frequency_change,absolute_shape_correlation,normalized_shape_l2_error\n";
    csv << std::setprecision(12);
    double max_freq_delta = 0.0, max_shape_error = 0.0, min_shape_corr = 1.0;
    for (size_t l = 1; l < levels.size(); ++l) {
        const auto& coarse = levels[l - 1];
        const auto& fine = levels[l];
        for (size_t m = 0; m < 10; ++m) {
            const double fc = coarse.modes.modes[m].frequency_hz;
            const double ff = fine.modes.modes[m].frequency_hz;
            const double freq_delta = std::abs(ff - fc) / std::abs(ff);
            const auto shape = compareShape(coarse.modes.modes[m], coarse.profile_nodes,
                                            fine.modes.modes[m], fine.profile_nodes);
            max_freq_delta = std::max(max_freq_delta, freq_delta);
            max_shape_error = std::max(max_shape_error, shape.error);
            min_shape_corr = std::min(min_shape_corr, shape.correlation);
            csv << coarse.profile_nodes << ',' << fine.profile_nodes << ',' << m + 1 << ','
                << fc << ',' << ff << ',' << freq_delta << ',' << shape.correlation << ','
                << shape.error << '\n';
        }
    }
    const auto summary_path = out_path.parent_path() / "mesh_convergence_summary.txt";
    std::ofstream summary(summary_path);
    if (!summary) { std::cerr << "Cannot write " << summary_path << '\n'; return 3; }
    summary << std::setprecision(12)
            << "geometry=corrected_inner_normal_offset\n"
            << "method=internal_axisymmetric_Q4_FEM\n"
            << "refinement=meridional_profile_only; one Q4 element across wall thickness\n"
            << "modes_compared=10\n"
            << "profile_nodes_sequence=41,51,61,71,81\n"
            << "fundamental_frequency_hz_by_level\n";
    for (const auto& level : levels)
        summary << level.profile_nodes << ',' << level.modes.modes.front().frequency_hz << '\n';
    const auto& last_coarse = levels[levels.size() - 2];
    const auto& last_fine = levels.back();
    double final_max_freq_delta = 0.0, final_max_shape_error = 0.0, final_min_shape_corr = 1.0;
    for (size_t m = 0; m < 10; ++m) {
        const double fc = last_coarse.modes.modes[m].frequency_hz;
        const double ff = last_fine.modes.modes[m].frequency_hz;
        final_max_freq_delta = std::max(final_max_freq_delta, std::abs(ff - fc) / std::abs(ff));
        const auto shape = compareShape(last_coarse.modes.modes[m], last_coarse.profile_nodes,
                                        last_fine.modes.modes[m], last_fine.profile_nodes);
        final_max_shape_error = std::max(final_max_shape_error, shape.error);
        final_min_shape_corr = std::min(final_min_shape_corr, shape.correlation);
    }
    summary << "\nfinal_adjacent_refinement=71_to_81\n"
            << "max_frequency_change_modes_1_to_10_percent=" << final_max_freq_delta * 100.0 << '\n'
            << "min_absolute_shape_correlation_modes_1_to_10=" << final_min_shape_corr << '\n'
            << "max_normalized_shape_l2_error_modes_1_to_10_percent=" << final_max_shape_error * 100.0 << '\n'
            << "all_level_pairs_max_frequency_change_percent=" << max_freq_delta * 100.0 << '\n'
            << "all_level_pairs_min_absolute_shape_correlation=" << min_shape_corr << '\n'
            << "all_level_pairs_max_normalized_shape_l2_error_percent=" << max_shape_error * 100.0 << '\n'
            << "interpretation=preliminary_convergence_only; refine_thickness_direction_before_final_mesh_independence_claim\n";
    std::cout << "Maximum adjacent-level frequency change (modes 1-10): "
              << max_freq_delta * 100.0 << "%\n"
              << "Minimum adjacent-level absolute shape correlation: " << min_shape_corr << '\n'
              << "Maximum adjacent-level normalized shape L2 error: " << max_shape_error * 100.0 << "%\n"
              << "Final 71->81 max frequency change: " << final_max_freq_delta * 100.0 << "%\n"
              << "Final 71->81 min shape correlation: " << final_min_shape_corr << '\n'
              << "Results: " << out_path.string() << '\n';
    return csv && summary ? 0 : 4;
}
