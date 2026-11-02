#include "core/bell_types.h"
#include "core/modal_types.hpp"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    try {
        BellGeometry geometry;
        geometry.H = 1000.0;
        geometry.R_mouth = 500.0;
        geometry.t_top = 10.0;
        geometry.t_bottom = 30.0;
        geometry.r_knots = {0.0, 0.13, 0.22, 0.46, 1.0};
        geometry.ellipticity = 1.0;
        BellMaterial material;
        BellMeshConfig meshConfig;
        meshConfig.nNodes = 15;
        const auto model = make_cylindrical_bell_model(geometry, material, meshConfig);
        require(model.mesh_si.inner_profile.size() == 15, "inner geometry profile missing");

        Excitation::MeridionalInfluence influence;
        CylindricalModalEvaluationConfig config;
        config.harmonics = {0, 1};
        config.angular_samples = 32;
        config.requested_modes = 8;
        config.max_iterations = 5000;
        config.tolerance = 1.0e-6;
        ModalData modes;
        require(evaluate_cylindrical_bell_modes(model, influence.make(model.mesh_si), config, modes),
                "direct FEM modal evaluation failed");

        const std::size_t outerCount = model.mesh_si.outer_profile.size();
        const std::size_t innerCount = model.mesh_si.inner_profile.size();
        require(outerCount + innerCount <= model.mesh_si.nodes.size(), "invalid mesh boundary node layout");
        for (const auto& mode : modes.modes) {
            require(mode.inner_surface_radial.size() == innerCount &&
                    mode.inner_surface_axial.size() == innerCount,
                    "scalar inner boundary export has wrong size");
            require(mode.inner_surface_harmonics == mode.outer_surface_harmonics,
                    "inner and outer Fourier bases differ");
            const std::size_t basisCount = mode.inner_surface_harmonics.size();
            require(mode.inner_surface_fourier_radial.size() == innerCount * basisCount &&
                    mode.inner_surface_fourier_azimuthal.size() == innerCount * basisCount &&
                    mode.inner_surface_fourier_axial.size() == innerCount * basisCount,
                    "inner Fourier export has wrong size");
            require(mode.fem_harmonics.size() == basisCount, "FEM harmonic labels have unexpected size");
            const std::size_t dpn = 3 * basisCount;
            for (std::size_t k = 0; k < innerCount; ++k) {
                const std::size_t node = outerCount + k;
                for (std::size_t b = 0; b < basisCount; ++b) {
                    const std::size_t dst = k * basisCount + b;
                    const std::size_t srcR = node * dpn + b;
                    const std::size_t srcP = node * dpn + basisCount + b;
                    const std::size_t srcZ = node * dpn + 2 * basisCount + b;
                    require(std::abs(mode.inner_surface_fourier_radial[dst] - mode.fem_displacement[srcR]) < 1e-12,
                            "inner radial boundary differs from full FEM vector");
                    require(std::abs(mode.inner_surface_fourier_azimuthal[dst] - mode.fem_displacement[srcP]) < 1e-12,
                            "inner azimuthal boundary differs from full FEM vector");
                    require(std::abs(mode.inner_surface_fourier_axial[dst] - mode.fem_displacement[srcZ]) < 1e-12,
                            "inner axial boundary differs from full FEM vector");
                }
                require(std::abs(mode.inner_surface_radial[k] - mode.inner_surface_fourier_radial[k * basisCount]) < 1e-12,
                        "legacy inner radial view is not the zero-harmonic coefficient");
                require(std::abs(mode.inner_surface_axial[k] - mode.inner_surface_fourier_axial[k * basisCount]) < 1e-12,
                        "legacy inner axial view is not the zero-harmonic coefficient");
            }
        }
        std::cout << "inner_surface_modal_export_contract: PASS modes=" << modes.modes.size() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "inner_surface_modal_export_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
