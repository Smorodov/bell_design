#include "fem/cylindrical_bell_model.hpp"
#include "excitation/chinese_bell_strike.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Metrics {
    double integrated_weight = 0.0;
    double centroid_z = 0.0;
};

Metrics evaluate(int nodes, double width) {
    BellGeometry geometry;
    geometry.H = 1000.0;
    geometry.R_mouth = 500.0;
    geometry.t_top = 10.0;
    geometry.t_bottom = 30.0;
    geometry.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    geometry.ellipticity = 0.85;

    BellMaterial material;
    BellMeshConfig mesh;
    mesh.nNodes = nodes;
    const auto model = BellOptimizer::FEM::make_cylindrical_bell_model(
        geometry, material, mesh);

    const std::vector<int> harmonics{0, 2, 4};
    BellOptimizer::Excitation::ChineseBellStrikeConfig strike;
    strike.center_z = 0.62;
    strike.azimuth_rad = 0.73;
    strike.meridional_width = width;
    strike.angular_width_rad = 0.20;
    const auto load = BellOptimizer::Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, harmonics, geometry.ellipticity, strike, 128);
    if (load.dof_values.empty()) throw std::runtime_error("empty strike load");

    constexpr std::size_t basis_count = 5;
    Metrics out;
    double weighted_z = 0.0;
    for (std::size_t i = 0; i < model.mesh_si.outer_profile.size(); ++i) {
        const std::size_t base = i * 3u * basis_count;
        const double nr = load.dof_values[base];
        const double nphi = load.dof_values[base + 1u];
        const double nz = load.dof_values[base + 2u];
        const double weight = std::sqrt(nr * nr + nphi * nphi + nz * nz);
        out.integrated_weight += weight;
        weighted_z += weight * model.mesh_si.outer_profile[i].z;
    }
    if (!(out.integrated_weight > 0.0)) throw std::runtime_error("zero integrated contact weight");
    out.centroid_z = weighted_z / out.integrated_weight;
    return out;
}
}

int main() {
    const std::vector<int> meshes{15, 23, 31, 39, 47};
    const std::vector<double> widths{0.015, 0.035, 0.070};
    std::cout << std::setprecision(10);
    for (double width : widths) {
        std::vector<Metrics> values;
        for (int n : meshes) values.push_back(evaluate(n, width));
        const Metrics& fine = values.back();
        double previous_error = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < values.size(); ++i) {
            const double norm_error = std::abs(values[i].integrated_weight - 1.0);
            const double centroid_error = std::abs(values[i].centroid_z - fine.centroid_z);
            if (norm_error > 1.0e-12) throw std::runtime_error("contact weights are not normalized");
            if (i > 0 && centroid_error > previous_error * 1.25 + 1.0e-12)
                throw std::runtime_error("contact centroid refinement is unstable");
            previous_error = centroid_error;
            std::cout << "width_m=" << width << " nodes=" << meshes[i]
                      << " weight_sum=" << values[i].integrated_weight
                      << " centroid_z_m=" << values[i].centroid_z
                      << " fine_centroid_error_m=" << centroid_error << '\n';
        }
    }
    std::cout << "chinese_bell_strike_mesh_convergence: PASS\n";
    return 0;
}
