#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    BellMaterial material;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) {
        std::cerr << "mesh generation failed\n";
        return 2;
    }
    BellOptimizer::InternalModalSolver solver(true, 10);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "eigenpair_integrity", &modal_influence_1};
    BellOptimizer::ModalData result;
    if (!solver.solve(problem, result) || result.modes.size() < 8) {
        std::cerr << "modal solve failed or returned too few modes\n";
        return 3;
    }

    double max_residual = 0.0;
    double max_mass_error = 0.0;
    double max_orthogonality = 0.0;
    std::cout << std::setprecision(12)
              << "mode,frequency_hz,relative_residual,modal_mass,mass_orthogonality_error\n";
    for (size_t i = 0; i < result.modes.size(); ++i) {
        const auto& mode = result.modes[i];
        if (!std::isfinite(mode.eigen_residual_relative) ||
            !std::isfinite(mode.modal_mass) ||
            !std::isfinite(mode.mass_orthogonality_error)) return 4;
        max_residual = std::max(max_residual, mode.eigen_residual_relative);
        max_mass_error = std::max(max_mass_error, std::abs(mode.modal_mass - 1.0));
        max_orthogonality = std::max(max_orthogonality, mode.mass_orthogonality_error);
        std::cout << i + 1 << ',' << mode.frequency_hz << ','
                  << mode.eigen_residual_relative << ',' << mode.modal_mass << ','
                  << mode.mass_orthogonality_error << '\n';
    }
    std::cout << "max_relative_eigen_residual=" << max_residual << '\n'
              << "max_modal_mass_error=" << max_mass_error << '\n'
              << "max_mass_orthogonality_error=" << max_orthogonality << '\n';

    if (max_residual > 1e-6) return 5;
    if (max_mass_error > 1e-8) return 6;
    if (max_orthogonality > 1e-7) return 7;
    return 0;
}
