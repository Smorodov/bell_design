#pragma once

#include <vector>

namespace BellOptimizer::Excitation {

// A single linear structural mode coupled to a unilateral point contact.
// mode_shape_at_contact maps modal coordinate to displacement at the strike point.
struct ContactModalMode {
    double frequency_hz = 0.0;
    double modal_mass = 0.0;
    double mode_shape_at_contact = 0.0;
};

struct HammerContactConfig {
    double striker_mass_kg = 0.0;
    double initial_velocity_m_per_s = 0.0; // positive toward the bell
    double hertz_stiffness_n_per_m_1_5 = 0.0;
    double integration_step_seconds = 0.0;
    double duration_seconds = 0.0;
};

struct HammerContactResult {
    double initial_energy_joules = 0.0;
    double final_energy_joules = 0.0;
    double maximum_relative_energy_error = 0.0;
    double contact_duration_seconds = 0.0;
    std::vector<double> striker_displacement_m;
    std::vector<double> striker_velocity_m_per_s;
    std::vector<double> contact_force_n;
    std::vector<double> total_energy_joules;
    // [mode][time sample], modal coordinates in the supplied normalization.
    std::vector<std::vector<double>> modal_displacement;
    std::vector<std::vector<double>> modal_velocity;
};

// Conservative Hertz point contact (F=k*delta^(3/2), delta>0), coupled to
// undamped linear modes and a freely moving striker. No acoustic radiation,
// contact damping, plasticity, or finite contact-patch model is implied.
bool simulateHertzHammerContact(const HammerContactConfig& config,
                                const std::vector<ContactModalMode>& modes,
                                HammerContactResult& result);

} // namespace BellOptimizer::Excitation
