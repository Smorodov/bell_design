#include "bell_generator.h"
#include "core/bell_types.h"
#include "core/modal_tracking.hpp"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace BellOptimizer;

namespace {
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> values;
    std::stringstream stream(s);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    return values;
}

BellGeometry parseGeometry(const std::vector<std::string>& fields) {
    BellGeometry geometry;
    geometry.H = std::stod(fields[2]);
    geometry.R_mouth = std::stod(fields[3]);
    geometry.t_top = std::stod(fields[4]);
    geometry.t_bottom = std::stod(fields[5]);
    geometry.r_knots = {0.0, std::stod(fields[6]), std::stod(fields[7]),
                        std::stod(fields[8]), 1.0};
    geometry.ellipticity = std::stod(fields[9]);
    return geometry;
}

bool sameMeshOrdering(const BellMeshData& legacy, const BellMeshData& direct) {
    if (legacy.nodes.size() != direct.nodes.size()) return false;
    for (std::size_t i = 0; i < legacy.nodes.size(); ++i) {
        // BellGenerator stores millimetres; the direct cylindrical model stores SI.
        if (std::abs(legacy.nodes[i].r * 1e-3 - direct.nodes[i].r) > 1e-10 ||
            std::abs(legacy.nodes[i].z * 1e-3 - direct.nodes[i].z) > 1e-10)
            return false;
    }
    return true;
}

double fullVolumeLumpedMassOverlap(const ModalMode& legacy, const ModalMode& direct) {
    const std::size_t nodes = legacy.nodal_radial.size();
    if (nodes == 0 || legacy.nodal_axial.size() != nodes ||
        direct.fem_displacement.size() != nodes * 3 ||
        direct.fem_mass_lumped_weights.size() != nodes * 3 ||
        direct.fem_harmonics != std::vector<int>{0}) return -1.0;

    double dot = 0.0, legacy_norm2 = 0.0, direct_norm2 = 0.0;
    for (std::size_t node = 0; node < nodes; ++node) {
        for (const std::size_t component : {std::size_t(0), std::size_t(2)}) {
            const std::size_t dof = node * 3 + component;
            const double a = component == 0 ? legacy.nodal_radial[node]
                                             : legacy.nodal_axial[node];
            const double b = direct.fem_displacement[dof];
            const double weight = direct.fem_mass_lumped_weights[dof];
            if (!(weight >= 0.0) || !std::isfinite(weight)) return -1.0;
            dot += weight * a * b;
            legacy_norm2 += weight * a * a;
            direct_norm2 += weight * b * b;
        }
    }
    if (!(legacy_norm2 > 0.0) || !(direct_norm2 > 0.0)) return -1.0;
    return std::abs(dot) / std::sqrt(legacy_norm2 * direct_norm2);
}

std::vector<int> maximumAssignment(const std::vector<std::vector<double>>& score) {
    const std::size_t rows = score.size();
    if (rows == 0) return {};
    const std::size_t cols = score.front().size();
    if (cols < rows || cols > 20) return {};
    for (const auto& row : score) if (row.size() != cols) return {};

    const std::size_t states = std::size_t(1) << cols;
    std::vector<double> dp(states, -std::numeric_limits<double>::infinity());
    std::vector<std::vector<int>> parent(rows, std::vector<int>(states, -1));
    dp[0] = 0.0;
    for (std::size_t row = 0; row < rows; ++row) {
        std::vector<double> next(states, -std::numeric_limits<double>::infinity());
        for (std::size_t mask = 0; mask < states; ++mask) {
            if (!std::isfinite(dp[mask])) continue;
            for (std::size_t col = 0; col < cols; ++col) {
                const std::size_t bit = std::size_t(1) << col;
                if (mask & bit) continue;
                const std::size_t next_mask = mask | bit;
                const double value = dp[mask] + score[row][col];
                if (value > next[next_mask]) {
                    next[next_mask] = value;
                    parent[row][next_mask] = static_cast<int>(col);
                }
            }
        }
        dp.swap(next);
    }

    std::size_t mask = 0;
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t state = 0; state < states; ++state) {
        if (std::isfinite(dp[state]) && dp[state] > best) {
            best = dp[state];
            mask = state;
        }
    }
    std::vector<int> assignment(rows, -1);
    for (std::size_t row = rows; row-- > 0;) {
        const int col = parent[row][mask];
        if (col < 0) return {};
        assignment[row] = col;
        mask &= ~(std::size_t(1) << static_cast<std::size_t>(col));
    }
    return assignment;
}

bool solveFast(const BellGeometry& geometry, const BellMaterial& material,
               const BellMeshConfig& mesh, ModalData& result, BellMeshData& mesh_data) {
    BellGenerator generator(material, mesh);
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    mesh_data = generator.getMeshData();
    Excitation::MeridionalInfluence influence;
    ModalProblem problem{mesh_data, material, ".", "round_full_volume_fast", &influence};
    InternalModalSolver solver(true, 8, true);
    return solver.solve(problem, result) && result.modes.size() >= 5;
}

bool solveDirect(const BellGeometry& geometry, const BellMaterial& material,
                 const BellMeshConfig& mesh, ModalData& result, BellMeshData& legacy_mesh) {
    const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh);
    if (!sameMeshOrdering(legacy_mesh, model.mesh_si)) return false;
    Excitation::MeridionalInfluence influence;
    FEM::CylindricalModalEvaluationConfig config;
    config.harmonics = {0};
    config.angular_samples = 32;
    config.requested_modes = 8;
    config.max_iterations = 1600;
    config.tolerance = 1e-6;
    config.minimum_frequency_hz = 1.0;
    return FEM::evaluate_cylindrical_bell_modes(
        model, influence.make(model.mesh_si), config, result) && result.modes.size() >= 5;
}
} // namespace

int main(int argc, char** argv) {
    try {
        const std::string input = argc > 1 ? argv[1] : "research/data/round_geometry_dataset.csv";
        const std::string output = argc > 2 ? argv[2] : "round_solver_full_volume_overlap.csv";
        const int requested_nodes = argc > 3 ? std::stoi(argv[3]) : 15;
        const std::size_t max_geometries = argc > 4 ? static_cast<std::size_t>(std::stoul(argv[4])) : 0;
        std::ifstream source(input);
        if (!source) throw std::runtime_error("cannot open geometry dataset: " + input);
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("cannot open output CSV: " + output);

        std::string line;
        std::getline(source, line);
        csv << "geometry,fast_mode,direct_mode,fast_frequency_hz,direct_frequency_hz,frequency_relative_error,lumped_mass_full_volume_overlap\n";
        csv << std::setprecision(12);

        BellMaterial material;
        BellMeshConfig mesh;
        mesh.nNodes = requested_nodes;
        if (mesh.nNodes < 5) throw std::runtime_error("mesh node count must be at least 5");
        std::size_t geometries = 0, failures = 0, pairs = 0;
        double minimum_overlap = 1.0, maximum_frequency_error = 0.0;

        while (std::getline(source, line)) {
            if (max_geometries > 0 && geometries >= max_geometries) break;
            if (line.empty()) continue;
            const auto fields = split(line);
            if (fields.size() < 17 || fields[16] != "PASS") continue;
            const BellGeometry geometry = parseGeometry(fields);
            ModalData fast, direct;
            BellMeshData legacy_mesh;
            if (!solveFast(geometry, material, mesh, fast, legacy_mesh) ||
                !solveDirect(geometry, material, mesh, direct, legacy_mesh)) {
                ++failures;
                std::cerr << "SOLVE_OR_MESH_FAIL " << fields[1] << '\n';
                continue;
            }

            const std::size_t fast_count = std::min<std::size_t>(5, fast.modes.size());
            const std::size_t direct_count = std::min<std::size_t>(8, direct.modes.size());
            std::vector<std::vector<double>> score(fast_count, std::vector<double>(direct_count, -1.0));
            for (std::size_t i = 0; i < fast_count; ++i)
                for (std::size_t j = 0; j < direct_count; ++j)
                    score[i][j] = fullVolumeLumpedMassOverlap(fast.modes[i], direct.modes[j]);
            const auto assignment = maximumAssignment(score);
            if (assignment.size() != fast_count) {
                ++failures;
                std::cerr << "ASSIGNMENT_FAIL " << fields[1] << '\n';
                continue;
            }

            ++geometries;
            for (std::size_t i = 0; i < fast_count; ++i) {
                const std::size_t j = static_cast<std::size_t>(assignment[i]);
                const double overlap = score[i][j];
                const double frequency_error = std::abs(fast.modes[i].frequency_hz -
                    direct.modes[j].frequency_hz) / std::max(direct.modes[j].frequency_hz, 1e-12);
                minimum_overlap = std::min(minimum_overlap, overlap);
                maximum_frequency_error = std::max(maximum_frequency_error, frequency_error);
                ++pairs;
                csv << fields[1] << ',' << i << ',' << j << ','
                    << fast.modes[i].frequency_hz << ',' << direct.modes[j].frequency_hz << ','
                    << frequency_error << ',' << overlap << '\n';
            }
            if (geometries % 10 == 0) std::cout << "completed=" << geometries << '\n';
        }

        std::cout << std::setprecision(12)
                  << "accepted_geometries=" << geometries << " failures=" << failures << '\n'
                  << "matched_first_five_pairs=" << pairs << '\n'
                  << "minimum_lumped_mass_full_volume_overlap=" << minimum_overlap << '\n'
                  << "maximum_matched_frequency_relative_error=" << maximum_frequency_error << '\n'
                  << "scope=m0_only; mesh_nodes=" << mesh.nNodes
                  << "; first-five-to-eight one-to-one rank-1 full-field assignment; diagonal row-sum mass weights\n";
        return geometries > 0 && failures == 0 && pairs == geometries * 5 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "round_solver_full_volume_probe: FAIL: " << error.what() << '\n';
        return 1;
    }
}
