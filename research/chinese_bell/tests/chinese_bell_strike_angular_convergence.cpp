#include "fem/cylindrical_bell_model.hpp"
#include "excitation/chinese_bell_strike.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace BellOptimizer;
    BellGeometry geometry;
    geometry.H = 1000.0;
    geometry.R_mouth = 500.0;
    geometry.t_top = 10.0;
    geometry.t_bottom = 30.0;
    geometry.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    geometry.ellipticity = 0.85;
    BellMaterial material;
    BellMeshConfig mesh;
    mesh.nNodes = 23;
    const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh);
    const std::vector<int> harmonics{0, 2, 4};
    Excitation::ChineseBellStrikeConfig strike;
    strike.center_z = 0.62;
    strike.azimuth_rad = 0.73;
    strike.meridional_width = 0.035;
    strike.angular_width_rad = 0.20;

    const auto l32 = Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, harmonics, geometry.ellipticity, strike, 32);
    const auto l64 = Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, harmonics, geometry.ellipticity, strike, 64);
    const auto l128 = Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, harmonics, geometry.ellipticity, strike, 128);
    if (l32.dof_values.size() != l64.dof_values.size() ||
        l64.dof_values.size() != l128.dof_values.size())
        throw std::runtime_error("angular quadrature changed load dimension");

    auto relative_difference = [](const auto& a, const auto& b) {
        double out = 0.0;
        for (std::size_t i = 0; i < a.dof_values.size(); ++i) {
            const double scale = std::max({1.0, std::abs(a.dof_values[i]), std::abs(b.dof_values[i])});
            out = std::max(out, std::abs(a.dof_values[i] - b.dof_values[i]) / scale);
        }
        return out;
    };
    const double d32_64 = relative_difference(l32, l64);
    const double d64_128 = relative_difference(l64, l128);
    if (!(d64_128 <= d32_64 * 1.10 + 1.0e-10))
        throw std::runtime_error("angular quadrature does not converge monotonically");
    if (!(d64_128 < 2.0e-3))
        throw std::runtime_error("angular quadrature error remains too large");

    std::cout << "chinese_bell_strike_angular_convergence: PASS\n"
              << "relative_change_32_to_64=" << d32_64 << '\n'
              << "relative_change_64_to_128=" << d64_128 << '\n';
    return 0;
}
