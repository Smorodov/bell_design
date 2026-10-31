#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_global_assembly.hpp"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::FEM;

    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 10;
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    require(generator.validateGeometry().passed && generator.buildMesh(), "production mesh generation failed");

    Excitation::MeridionalInfluence influence;
    ModalProblem problem{generator.getMeshData(), material, ".", "circular-limit-comparison", &influence};
    InternalModalSolver production(true, 8, true);
    ModalData reference;
    require(production.solve(problem, reference), "InternalModalSolver failed on reference mesh");
    require(reference.modes.size() >= 5, "InternalModalSolver returned too few modes");

    // Production units are mm, MPa and tonne/mm^3. Convert exactly to SI for
    // the new vector formulation: metres, pascals and kg/m^3.
    BellMeshData si_mesh = generator.getMeshData();
    for (auto& node : si_mesh.nodes) { node.r *= 1.0e-3; node.z *= 1.0e-3; }
    const double E_pa = material.E * 1.0e6;
    const double rho_kg_m3 = material.rho * 1.0e12;
    const auto system = assemble_cylindrical_global_system(si_mesh, E_pa, material.nu,
                                                            rho_kg_m3, {0}, 32);

    // Match the production meridional subspace: constrain radial displacement
    // on the axis, and suppress the decoupled n=0 circumferential torsion DOF.
    // Axial displacement remains free at the axis in both formulations.
    std::vector<int> fixed;
    for (std::size_t i = 0; i < si_mesh.nodes.size(); ++i) {
        fixed.push_back(static_cast<int>(i * 3u + 1u)); // uphi = 0 for all nodes
        if (std::abs(si_mesh.nodes[i].r) < 1.0e-12)
            fixed.push_back(static_cast<int>(i * 3u)); // ur = 0 on the axis
    }
    const auto vector_result = solve_cylindrical_global_modes(system, fixed, 9, 1800, 5.0e-5, true);
    // The free body contains a zero-frequency rigid axial translation. The
    // current API marks the solve incomplete when this zero mode is filtered
    // from the requested positive-mode count; accept the returned positive
    // modes only when their residual remains within the explicit audit bound.
    require(vector_result.frequencies_hz.size() >= 5 &&
            std::isfinite(vector_result.relative_residual) && vector_result.relative_residual < 1.0e-4,
            "new vector FEM returned too few positive modes or excessive residual");

    std::cout << std::setprecision(10);
    std::cout << "mode,production_hz,vector_n0_hz,relative_error\n";
    double max_relative_error = 0.0;
    const std::size_t count = std::min<std::size_t>(8, std::min(reference.modes.size(), vector_result.frequencies_hz.size()));
    for (std::size_t i = 0; i < count; ++i) {
        const double a = reference.modes[i].frequency_hz;
        const double b = vector_result.frequencies_hz[i];
        const double error = std::abs(a - b) / std::max(std::abs(a), 1.0e-12);
        max_relative_error = std::max(max_relative_error, error);
        std::cout << (i + 1) << ',' << a << ',' << b << ',' << error << '\n';
    }
    std::cout << "max_relative_error=" << max_relative_error
              << " vector_residual=" << vector_result.relative_residual << '\n';

    // This is deliberately a diagnostic gate, not a hard claim of production
    // equivalence: the two solvers use different eigensolver paths and meshes
    // are shared but material/coordinate conversion and axis handling are explicit.
    require(max_relative_error < 0.08,
            "circular-limit frequencies diverge from the production meridional FEM by more than 8%");
    return 0;
}
