#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <cmath>
#include <iostream>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;

    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 2;

    BellOptimizer::InternalModalSolver solver(true, 10);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "acoustic", &modal_influence_1};
    if (!solver.solve(problem, data) || data.modes.size() < 5) return 3;

    size_t nonzero = 0;
    double max_abs_qd = 0.0;
    for (const auto& mode : data.modes) {
        if (!std::isfinite(mode.acoustic_volume_displacement)) return 4;
        if (!std::isfinite(mode.modal_excitation)) return 5;
        max_abs_qd = std::max(max_abs_qd, std::abs(mode.acoustic_volume_displacement));
        if (std::abs(mode.acoustic_volume_displacement) > 1e-14) ++nonzero;
    }

    std::cout << "modes=" << data.modes.size() << "\n";
    std::cout << "nonzero_acoustic_coefficients=" << nonzero << "\n";
    std::cout << "max_abs_acoustic_volume_displacement=" << max_abs_qd << "\n";

    // The surface-radiation coefficient is a mode-shape property. It must not
    // be silently replaced by effective modal mass when the shape is present.
    return nonzero >= 3 ? 0 : 6;
}
