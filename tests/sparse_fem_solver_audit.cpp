#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <cmath>
#include <iostream>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 75; // 298 reduced DOFs: forces the sparse iterative production path.
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 1;

    BellOptimizer::InternalModalSolver solver(true, 5, true);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "sparse_fem_solver_audit", &modal_influence_1};
    BellOptimizer::ModalData result;
    if (!solver.solve(problem, result)) return 2;
    if (result.modes.size() != 5 || !result.matrix_audit.available || result.matrix_audit.reduced_dof <= 128) return 3;
    if (result.matrix_audit.assembly_dense_bytes != 0) return 4;
    for (const auto& mode : result.modes) {
        if (!(mode.frequency_hz > 0.0) || !std::isfinite(mode.frequency_hz)) return 5;
        if (!(mode.eigen_residual_relative < 1e-5)) return 6;
        if (std::abs(mode.modal_mass - 1.0) > 1e-8) return 7;
        if (!std::isfinite(mode.mass_orthogonality_error) || mode.mass_orthogonality_error > 1e-4) return 8;
    }
    std::cout << "sparse_fem_solver_audit=PASS\n"
              << "reduced_dof=" << result.matrix_audit.reduced_dof << '\n'
              << "K_nnz=" << result.matrix_audit.stiffness.exact_nonzeros << '\n'
              << "M_nnz=" << result.matrix_audit.mass.exact_nonzeros << '\n'
              << "first_frequency_hz=" << result.modes.front().frequency_hz << '\n'
              << "max_residual=";
    double max_residual = 0.0;
    for (const auto& mode : result.modes) max_residual = std::max(max_residual, mode.eigen_residual_relative);
    std::cout << max_residual << '\n';
    return 0;
}
