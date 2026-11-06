#pragma once

#include "clapper_dynamics.hpp"
#include <vector>

namespace BellOptimizer::Excitation {

struct ClapperModalEvaluationConfig {
    ClapperDynamicsConfig clapper;
    // Explicit axial reference points in the same SI coordinate system as profile.
    // The fraction is zero at top_z and one at mouth_z, independent of axis sign.
    double top_z_m = 0.0;
    double mouth_z_m = 1.0;
    double admissible_strike_fraction_min = 0.50;
    double admissible_strike_fraction_max = 0.70;
};

struct ClapperModalEvaluationResult {
    ClapperDynamicsResult contact;
    double first_contact_z_m = 0.0;
    double strike_height_fraction = 0.0;
    bool strike_height_admissible = false;
    // Peak free-response pressure proxy per mode, preserving each mode's
    // signed acoustic coefficient only through the phase pair below.
    std::vector<double> modal_pressure_cosine;
    std::vector<double> modal_pressure_sine;
    std::vector<double> modal_pressure_magnitude;
};

// Evaluates one physically specified clapper strike against a self-consistent
// modal state. Input profile and all lengths must be SI. No implicit unit conversion.
bool evaluateClapperModalState(
    const ClapperModalEvaluationConfig& config,
    const std::vector<Point2D>& inner_profile_si,
    const std::vector<ModalMode>& modes,
    ClapperModalEvaluationResult& result);

} // namespace BellOptimizer::Excitation
