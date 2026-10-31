#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 61;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 1;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "sparse_dense_mode_comparison", &modal_influence_1};
    BellOptimizer::ModalData dense, sparse;
    BellOptimizer::InternalModalSolver dense_solver(true, 5, true, 1000000);
    BellOptimizer::InternalModalSolver sparse_solver(true, 5, true, 0);
    if (!dense_solver.solve(problem, dense) || !sparse_solver.solve(problem, sparse)) return 2;
    if (dense.modes.size() != sparse.modes.size()) return 3;
    double max_frequency_error = 0.0;
    double min_shape_correlation = 1.0;
    double max_sparse_residual = 0.0;
    for (std::size_t i = 0; i < dense.modes.size(); ++i) {
        const auto& a = dense.modes[i];
        const auto& b = sparse.modes[i];
        const double frequency_error = std::abs(a.frequency_hz - b.frequency_hz) /
            std::max(std::abs(a.frequency_hz), 1e-30);
        max_frequency_error = std::max(max_frequency_error, frequency_error);
        if (a.nodal_radial.size() != b.nodal_radial.size() || a.nodal_axial.size() != b.nodal_axial.size()) return 4;
        double dot = 0.0, norm_a = 0.0, norm_b = 0.0;
        for (std::size_t j = 0; j < a.nodal_radial.size(); ++j) {
            dot += a.nodal_radial[j] * b.nodal_radial[j] + a.nodal_axial[j] * b.nodal_axial[j];
            norm_a += a.nodal_radial[j] * a.nodal_radial[j] + a.nodal_axial[j] * a.nodal_axial[j];
            norm_b += b.nodal_radial[j] * b.nodal_radial[j] + b.nodal_axial[j] * b.nodal_axial[j];
        }
        const double correlation = std::abs(dot) / std::max(std::sqrt(norm_a * norm_b), 1e-30);
        min_shape_correlation = std::min(min_shape_correlation, correlation);
        max_sparse_residual = std::max(max_sparse_residual, b.eigen_residual_relative);
    }
    if (max_frequency_error > 2e-4 || min_shape_correlation < 0.995 || max_sparse_residual > 1e-5) return 5;
    std::cout << "sparse_dense_mode_comparison=PASS\n"
              << "profile_nodes=" << config.nNodes << " reduced_dof=" << dense.matrix_audit.reduced_dof << '\n'
              << "max_frequency_relative_error=" << max_frequency_error << '\n'
              << "min_euclidean_shape_correlation=" << min_shape_correlation << '\n'
              << "max_sparse_eigen_residual=" << max_sparse_residual << '\n';
    return 0;
}
