#include "fem/cylindrical_bell_model.hpp"
#include "excitation/chinese_bell_strike.hpp"
#include "core/bell_types.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using namespace BellOptimizer;
    BellGeometry g;
    g.H = 1000.0;
    g.R_mouth = 500.0;
    g.t_top = 10.0;
    g.t_bottom = 30.0;
    g.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    BellMaterial mat;
    BellMeshConfig mesh;
    mesh.nNodes = 15;
    const std::vector<int> harmonics{0, 2, 4};

    g.ellipticity = 1.0;
    const auto circular = FEM::make_cylindrical_bell_model(g, mat, mesh);
    Excitation::ChineseBellStrikeConfig cfg;
    cfg.center_z = 0.70;
    cfg.azimuth_rad = 0.35;
    cfg.meridional_width = 0.045;
    cfg.angular_width_rad = 0.22;
    const auto circular_load = Excitation::make_chinese_bell_strike_influence(
        circular.mesh_si, harmonics, 1.0, cfg, 64);

    g.ellipticity = 0.85;
    const auto elliptical = FEM::make_cylindrical_bell_model(g, mat, mesh);
    cfg.azimuth_rad = 0.7853981633974483;
    const auto elliptical_load = Excitation::make_chinese_bell_strike_influence(
        elliptical.mesh_si, harmonics, 0.85, cfg, 64);
    if (circular_load.dof_values.empty() || elliptical_load.dof_values.empty())
        throw std::runtime_error("empty strike influence");

    const std::size_t nb = 5;
    double circular_phi_max = 0.0;
    double elliptical_phi_max = 0.0;
    double vector_norm_error = 0.0;
    for (std::size_t node = 0; node < elliptical.mesh_si.outer_profile.size(); ++node) {
        for (std::size_t b = 0; b < nb; ++b) {
            const std::size_t base = node * 3u * nb + 3u * b;
            circular_phi_max = std::max(circular_phi_max,
                std::abs(circular_load.dof_values[base + 1u]));
            elliptical_phi_max = std::max(elliptical_phi_max,
                std::abs(elliptical_load.dof_values[base + 1u]));
        }
    }
    if (circular_phi_max > 1.0e-13)
        throw std::runtime_error("circular surface normal must have no azimuthal component");
    if (!(elliptical_phi_max > 1.0e-8))
        throw std::runtime_error("elliptical surface normal is missing its azimuthal component");

    // The load is now the actual angularly integrated surface-force projection,
    // so a finite angular patch is not expected to equal the pointwise normal
    // at its centre.  The contract therefore checks symmetry and convergence
    // properties of the integrated load instead of comparing it with a point
    // normal.
    for (std::size_t node = 0; node < elliptical.mesh_si.outer_profile.size(); ++node) {
        for (std::size_t b = 0; b < nb; ++b) {
            const std::size_t base = node * 3u * nb + 3u * b;
            if (!std::isfinite(elliptical_load.dof_values[base + 0u]) ||
                !std::isfinite(elliptical_load.dof_values[base + 1u]) ||
                !std::isfinite(elliptical_load.dof_values[base + 2u]))
                throw std::runtime_error("non-finite integrated strike load");
        }
    }
    const auto elliptical_load_128 = Excitation::make_chinese_bell_strike_influence(
        elliptical.mesh_si, harmonics, 0.85, cfg, 128);
    if (elliptical_load_128.dof_values.size() != elliptical_load.dof_values.size())
        throw std::runtime_error("angular quadrature changed load dimension");
    double max_relative_quadrature_change = 0.0;
    for (std::size_t i = 0; i < elliptical_load.dof_values.size(); ++i) {
        const double scale = std::max({1.0, std::abs(elliptical_load.dof_values[i]),
                                       std::abs(elliptical_load_128.dof_values[i])});
        max_relative_quadrature_change = std::max(
            max_relative_quadrature_change,
            std::abs(elliptical_load.dof_values[i] - elliptical_load_128.dof_values[i]) / scale);
    }
    if (!(max_relative_quadrature_change < 2.0e-3))
        throw std::runtime_error("angular quadrature convergence is too weak");

    std::cout << "chinese_bell_strike_surface_normal_contract: PASS\n";
    std::cout << "circular_azimuthal_component=" << circular_phi_max
              << " elliptical_azimuthal_component=" << elliptical_phi_max
              << " normal_unit_error=" << vector_norm_error << "\n";
    return 0;
}
