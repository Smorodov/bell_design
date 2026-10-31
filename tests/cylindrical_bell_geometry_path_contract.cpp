#include "fem/cylindrical_bell_model.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static BellGeometry production_geometry(double q) {
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

        const auto circular = make_cylindrical_bell_model(production_geometry(1.0), material, config);
        const auto elliptical = make_cylindrical_bell_model(production_geometry(0.85), material, config);

        require(!circular.mesh_si.nodes.empty() && circular.mesh_si.nodes.size() == 30,
                "production BellGeometry path did not build the expected production mesh");
        require(circular.mesh_si.elements.size() == 14,
                "production BellGeometry path returned an unexpected element count");
        require(circular.young_modulus_pa == 1.1e10,
                "production Young modulus was not converted from MPa to Pa");
        require(std::abs(circular.density_kg_m3 - 8600.0) < 1e-9,
                "production density was not converted from tonne/mm^3 to kg/m^3");
        require(circular.mesh_si.nodes.back().z > 0.1 && circular.mesh_si.nodes.back().r > 0.01,
                "production mesh was not converted from millimetres to metres");
        require(std::abs(circular.ellipticity - 1.0) < 1e-15,
                "circular q was not preserved");
        require(std::abs(elliptical.ellipticity - 0.85) < 1e-15,
                "elliptical q was not preserved");

        // The adapter must be a pure integration layer: q=1 must produce exactly
        // the same matrices as the existing global FEM assembler when given the
        // same SI mesh/material data.
        const auto circle_adapter = assemble_cylindrical_bell_system(circular, {0, 1, 2, 3}, 64);
        const auto circle_direct = assemble_cylindrical_global_system(
            circular.mesh_si, circular.young_modulus_pa, circular.poisson_ratio,
            circular.density_kg_m3, {0, 1, 2, 3}, 64, 0.0, 1.0);
        require(circle_adapter.stiffness.n == circle_direct.stiffness.n,
                "adapter and direct circular stiffness dimensions differ");
        require(circle_adapter.mass.n == circle_direct.mass.n,
                "adapter and direct circular mass dimensions differ");
        require(circle_adapter.stiffness.values.size() == circle_direct.stiffness.values.size(),
                "adapter and direct circular stiffness sparsity differs");

        double max_adapter_delta = 0.0;
        for (std::size_t i = 0; i < circle_adapter.stiffness.values.size(); ++i)
            max_adapter_delta = std::max(max_adapter_delta,
                std::abs(circle_adapter.stiffness.values[i] - circle_direct.stiffness.values[i]));
        require(max_adapter_delta < 1e-12,
                "adapter changed the circular FEM stiffness matrix");

        // q must now enter the actual FEM assembly, not a post-processing
        // coordinate morphing step. Compare the real assembled matrices.
        const auto ellipse_system = assemble_cylindrical_bell_system(
            elliptical, {0, 1, 2, 3}, 64);
        require(ellipse_system.stiffness.n == circle_adapter.stiffness.n,
                "elliptical and circular FEM systems have different dimensions");
        double circle_k_l1 = 0.0, ellipse_k_l1 = 0.0;
        double circle_m_l1 = 0.0, ellipse_m_l1 = 0.0;
        for (double v : circle_adapter.stiffness.values) circle_k_l1 += std::abs(v);
        for (double v : ellipse_system.stiffness.values) ellipse_k_l1 += std::abs(v);
        for (double v : circle_adapter.mass.values) circle_m_l1 += std::abs(v);
        for (double v : ellipse_system.mass.values) ellipse_m_l1 += std::abs(v);
        const double relative_k_change = std::abs(ellipse_k_l1 - circle_k_l1) / std::max(circle_k_l1, 1e-300);
        const double relative_m_change = std::abs(ellipse_m_l1 - circle_m_l1) / std::max(circle_m_l1, 1e-300);
        require(relative_k_change > 1e-4, "q=0.85 did not materially change the assembled stiffness matrix");
        require(relative_m_change > 1e-4, "q=0.85 did not materially change the assembled mass matrix");

        std::cout << "cylindrical_bell_geometry_path_contract: PASS\n";
        std::cout << "mesh_nodes=" << circular.mesh_si.nodes.size()
                  << " elements=" << circular.mesh_si.elements.size()
                  << " dofs=" << circle_adapter.stiffness.n << "\n";
        std::cout << "E_pa=" << circular.young_modulus_pa
                  << " rho_kg_m3=" << circular.density_kg_m3
                  << " max_adapter_matrix_delta=" << max_adapter_delta << "\n";
        std::cout << "q085_relative_stiffness_change=" << relative_k_change
                  << " q085_relative_stiffness_change=" << relative_k_change
                  << " q085_relative_mass_change=" << relative_m_change << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "cylindrical_bell_geometry_path_contract: FAIL: " << e.what() << "\n";
        return 1;
    }
}
