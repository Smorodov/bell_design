#include "bell_generator.h"
#include "core/bell_types.h"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/internal_modal_solver.hpp"
#include "math/modal_coordinate_morphing.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static BellGeometry reference_geometry(double q) {
    BellGeometry g;
    g.H = 658.3175709034034;
    g.R_mouth = 355.9758748202059;
    g.t_top = 22.85675347179014;
    g.t_bottom = 46.62524667849419;
    g.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    g.ellipticity = q;
    return g;
}

int main() {
    try {
        BellMaterial material;
        BellMeshConfig config;
        config.nNodes = 15;

        auto geometry = reference_geometry(0.85);
        geometry.scale(1000.0 / geometry.H);
        const auto model = make_cylindrical_bell_model(geometry, material, config);
        Excitation::MeridionalInfluence influence;

        CylindricalModalEvaluationConfig settings;
        settings.harmonics = {0, 2, 4};
        settings.angular_samples = 64;
        settings.requested_modes = 12;
        settings.max_iterations = 5000;
        settings.tolerance = 1.0e-6;
        settings.minimum_frequency_hz = 1.0;

        ModalData direct;
        require(evaluate_cylindrical_bell_modes(
                    model, influence.make(model.mesh_si), settings, direct),
                "direct cylindrical FEM evaluation failed");
        require(direct.modes.size() >= 5, "direct FEM returned too few positive modes");
        for (const auto& mode : direct.modes) {
            require(std::isfinite(mode.frequency_hz) && mode.frequency_hz > 0.0,
                    "direct FEM returned invalid frequency");
            require(std::isfinite(mode.amplitude) && mode.amplitude >= 0.0,
                    "direct FEM returned invalid amplitude");
        }

        // The direct path must not be a hidden frequency correction. Its first
        // modes are taken from the actual q=0.85 K/M eigenproblem.
        BellGenerator generator(material, config);
        generator.setGeometry(geometry);
        require(generator.validateGeometry().passed && generator.buildMesh(),
                "reference production mesh failed");
        ModalProblem problem{generator.getMeshData(), material, ".", "direct_bridge_reference", &influence};
        InternalModalSolver legacy(true, 8, true);
        ModalData base;
        require(legacy.solve(problem, base), "legacy production solve failed");
        require(base.modes.size() >= 5, "legacy production solve returned too few modes");

        std::cout << "direct_modes=";
        for (std::size_t i = 0; i < std::min<std::size_t>(5, direct.modes.size()); ++i)
            std::cout << ' ' << direct.modes[i].frequency_hz;
        std::cout << "\n";
        std::cout << "legacy_circular_modes=";
        for (std::size_t i = 0; i < std::min<std::size_t>(5, base.modes.size()); ++i)
            std::cout << ' ' << base.modes[i].frequency_hz;
        std::cout << "\n";
        std::cout << "direct_elliptical_fem_evaluation: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "direct_elliptical_fem_evaluation: FAIL: " << e.what() << "\n";
        return 1;
    }
}
