#include <algorithm>
#include "excitation/clapper_dynamics.hpp"
#include "core/bell_material.h"

#include <cmath>
#include <iostream>

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::Excitation;
    ClapperScalingParams scaling;
    BellMaterial material = default_material_preset().material;
    ClapperDynamicsConfig one, small;
    ClapperScalingResult summary;
    if (!makeScaledClapperDynamicsConfig(1.0, 0.0, scaling, material, one, nullptr)) return 1;
    if (!makeScaledClapperDynamicsConfig(0.2, 0.0, scaling, material, small, &summary)) return 2;
    auto close = [](double a, double b, double rel=1e-10) {
        return std::abs(a-b) <= rel * std::max({1.0, std::abs(a), std::abs(b)});
    };
    if (!close(small.mass_kg, 0.04) || !close(small.centre_of_mass_length_m, 0.092) ||
        !close(small.tip_length_m, 0.104) || !close(small.tip_radius_m, 0.007) ||
        !close(small.moment_of_inertia_kg_m2, scaling.reference_moment_of_inertia_kg_m2 * std::pow(0.2, 5.0)) ||
        !close(summary.linear_scale, 0.2) ||
        !close(small.initial_angular_velocity_rad_s, 3.2 / std::sqrt(0.2))) {
        std::cerr << "clapper similarity scaling law mismatch\n";
        return 3;
    }
    if (!(small.hertz_stiffness_n_per_m_1_5 > 0.0) ||
        !close(small.hertz_stiffness_n_per_m_1_5 / one.hertz_stiffness_n_per_m_1_5, std::sqrt(0.2)) ||
        !close(one.mass_kg, 5.0) || !close(one.tip_radius_m, 0.035)) return 4;
    if (makeScaledClapperDynamicsConfig(0.01, 0.0, scaling, material, small, nullptr)) {
        std::cerr << "unreasonable scale must be rejected\n";
        return 5;
    }
    std::cout << "clapper_scaling_contract: PASS\n";
    return 0;
}
