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

    // The normal is unit length before scalar spatial/Fourier weighting. Verify
    // the analytic cylindrical normal formula at a representative profile node.
    const std::size_t i = elliptical.mesh_si.outer_profile.size() / 2;
    const auto& p = elliptical.mesh_si.outer_profile;
    const std::size_t j0 = i - 1, j1 = i + 1;
    const double dr = p[j1].r - p[j0].r;
    const double dz = p[j1].z - p[j0].z;
    const double phi = cfg.azimuth_rad;
    const double q = 0.85;
    const double invq2 = 1.0 / (q * q);
    const double d = std::cos(phi) * std::cos(phi) + std::sin(phi) * std::sin(phi) * invq2;
    const double alpha = 1.0 / std::sqrt(d);
    const double gphi = -std::sin(phi) * std::cos(phi) * (invq2 - 1.0) / d;
    const double nr = dz, nphi = -dz * gphi, nz = -alpha * dr;
    const double expected_norm = std::sqrt(nr * nr + nphi * nphi + nz * nz);
    const double expected_nr = nr / expected_norm;
    const double expected_nphi = nphi / expected_norm;
    const double expected_nz = nz / expected_norm;
    const std::size_t base = i * 3u * nb;
    const double actual_nr = elliptical_load.dof_values[base + 0u];
    const double actual_nphi = elliptical_load.dof_values[base + 1u];
    const double actual_nz = elliptical_load.dof_values[base + 2u];
    const double actual_scale = std::sqrt(actual_nr * actual_nr + actual_nphi * actual_nphi + actual_nz * actual_nz);
    if (!(actual_scale > 0.0)) throw std::runtime_error("missing representative normal load");
    vector_norm_error = std::max({
        std::abs(actual_nr / actual_scale - expected_nr),
        std::abs(actual_nphi / actual_scale - expected_nphi),
        std::abs(actual_nz / actual_scale - expected_nz)});
    if (vector_norm_error > 1.0e-12)
        throw std::runtime_error("computed load direction differs from geometric surface normal");

    std::cout << "chinese_bell_strike_surface_normal_contract: PASS\n";
    std::cout << "circular_azimuthal_component=" << circular_phi_max
              << " elliptical_azimuthal_component=" << elliptical_phi_max
              << " normal_unit_error=" << vector_norm_error << "\n";
    return 0;
}
