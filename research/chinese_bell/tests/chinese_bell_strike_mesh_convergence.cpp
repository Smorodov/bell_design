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
    double coefficient_norm = 0.0;
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
    double weight_sum = 0.0;
    double coefficient_norm2 = 0.0;
    for (std::size_t i = 0; i < model.mesh_si.outer_profile.size(); ++i) {
        const std::size_t base = i * 3u * basis_count;
        double node_norm2 = 0.0;
        for (std::size_t b = 0; b < basis_count; ++b) {
            for (std::size_t c = 0; c < 3u; ++c) {
                const double v = load.dof_values[base + 3u * b + c];
                if (!std::isfinite(v)) throw std::runtime_error("non-finite strike coefficient");
                node_norm2 += v * v;
            }
        }
        const double weight = std::sqrt(node_norm2);
        coefficient_norm2 += node_norm2;
        weight_sum += weight;
        weighted_z += weight * model.mesh_si.outer_profile[i].z;
    }
    out.coefficient_norm = std::sqrt(coefficient_norm2);
    if (!(out.coefficient_norm > 0.0) || !(weight_sum > 0.0))
        throw std::runtime_error("zero integrated contact coefficient norm");
    out.centroid_z = weighted_z / weight_sum;
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
            const double centroid_error = std::abs(values[i].centroid_z - fine.centroid_z);
            if (i > 0 && centroid_error > previous_error * 1.25 + 1.0e-12)
                throw std::runtime_error("contact centroid refinement is unstable");
            previous_error = centroid_error;
            std::cout << "width_m=" << width << " nodes=" << meshes[i]
                      << "coefficient_norm=" << values[i].coefficient_norm
                      << " centroid_z_m=" << values[i].centroid_z
                      << " fine_centroid_error_m=" << centroid_error << '\n';
        }
    }
    std::cout << "chinese_bell_strike_mesh_convergence: PASS\n";
    return 0;
}
