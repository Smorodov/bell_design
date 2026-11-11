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

struct ClapperScalingResult {
    double linear_scale = 1.0;
    double mass_kg = 0.0;
    double moment_of_inertia_kg_m2 = 0.0;
    double centre_of_mass_length_m = 0.0;
    double tip_length_m = 0.0;
    double tip_radius_m = 0.0;
    double hertz_stiffness_n_per_m_1_5 = 0.0;
};

// Derive a geometrically similar clapper from a reference bell/clapper design.
// Lengths scale as s, mass as s^3, inertia as s^5, gravity-driven angular
// velocity as s^-1/2, and Hertz stiffness as sqrt(s) through tip radius.
bool makeScaledClapperDynamicsConfig(
    double final_bell_height_m,
    double inner_top_z_m,
    const ClapperScalingParams& scaling,
    const BellMaterial& bell_material,
    ClapperDynamicsConfig& config,
    ClapperScalingResult* summary = nullptr);

struct ClapperDynamicsResult {
    bool contact_occurred = false;
    double contact_duration_seconds = 0.0;
    double first_contact_time_seconds = 0.0;
    double separation_time_seconds = 0.0;
    double peak_contact_force_n = 0.0;
    // First-contact point on the undeformed inner meridian (SI coordinates).
    double first_contact_surface_r_m = 0.0;
    double first_contact_surface_z_m = 0.0;
    double initial_energy_joules = 0.0;
    double final_energy_joules = 0.0; // total conserved energy at separation
    double clapper_energy_at_separation_joules = 0.0;
    double bell_modal_energy_joules = 0.0;
    double bell_energy_transfer_fraction = 0.0; // bell modal energy / initial clapper energy
    std::vector<double> modal_energy_at_separation_joules;
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
