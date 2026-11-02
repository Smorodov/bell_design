#pragma once

#include "core/modal_types.hpp"
#include "core/bell_types.h"
#include <vector>

namespace BellOptimizer::Excitation {

struct ClapperDynamicsConfig {
    double mass_kg = 5.0;
    double moment_of_inertia_kg_m2 = 1.058;
    double centre_of_mass_length_m = 0.46;
    double tip_length_m = 0.52;
    double tip_radius_m = 0.035;
    double pivot_r_m = 0.0;
    double pivot_z_m = 0.0;
    double initial_angle_rad = 0.0;
    double initial_angular_velocity_rad_s = 3.2;
    double hertz_stiffness_n_per_m_1_5 = 0.0;
    double strike_azimuth_rad = 0.0;
    double integration_step_seconds = 2.0e-6;
    double maximum_duration_seconds = 0.25;
};

struct ClapperDynamicsResult {
    bool contact_occurred = false;
    double contact_duration_seconds = 0.0;
    double first_contact_time_seconds = 0.0;
    double separation_time_seconds = 0.0;
    double peak_contact_force_n = 0.0;
    double initial_energy_joules = 0.0;
    double final_energy_joules = 0.0;
    double maximum_relative_energy_error = 0.0;
    double separation_angle_rad = 0.0;
    double separation_angular_velocity_rad_s = 0.0;
    std::vector<double> modal_displacement_at_separation;
    std::vector<double> modal_velocity_at_separation;
    std::vector<double> contact_force_history_n;
    std::vector<double> angle_history_rad;
};

// Conservative planar physical-pendulum model coupled to the bell's inner-wall
// Fourier modal fields through a spherical Hertz contact. The clapper's mass,
// COM distance and pivot inertia define its rigid-body dynamics. The bell modes
// are undamped during the short impact; acoustic decay is applied after release.
bool simulateInternalClapperContact(
    const ClapperDynamicsConfig& config,
    const std::vector<Point2D>& inner_profile_si,
    const std::vector<ModalMode>& modes,
    ClapperDynamicsResult& result);

} // namespace BellOptimizer::Excitation
