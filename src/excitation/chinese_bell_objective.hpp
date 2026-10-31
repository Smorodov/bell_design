#pragma once

#include "core/modal_types.hpp"

namespace BellOptimizer::Excitation {

struct ChineseRadiationConfig {
    double air_density_kg_m3 = 1.20;
    double sound_speed_m_s = 343.0;
    double observation_distance_m = 1.0;
    // Direction cosine relative to the bell axis, in [-1,1].
    double axial_direction_cosine = 0.0;
};

struct ChineseRadiationResult {
    double monopole_m3 = 0.0;
    double dipole_m3 = 0.0;
    double pressure_pa_per_modal_coordinate = 0.0;
    double t60_seconds = 0.0;
    bool valid = false;
};

struct ChineseBellObjectiveConfig {
    double target_interval_cents = 1200.0;
    double interval_tolerance_cents = 35.0;
    double minimum_t60_seconds = 0.5;
    double maximum_competitor_ratio = 0.35;
    ChineseRadiationConfig radiation;
};

struct ChineseBellObjectiveResult {
    double interval_cents = 0.0;
    double interval_error_cents = 0.0;
    double long_strike_radiated_pressure = 0.0;
    double short_strike_radiated_pressure = 0.0;
    double long_branch_selectivity = 0.0;
    double short_branch_selectivity = 0.0;
    double weaker_target_t60_seconds = 0.0;
    bool interval_gate = false;
    bool excitation_gate = false;
    bool radiation_gate = false;
    bool decay_gate = false;
    bool passed = false;
};

// Far-field monopole + first axial-dipole proxy for a single mass-normalized mode.
// This is a directional linear-acoustics estimate, not a calibrated SPL prediction.
ChineseRadiationResult evaluate_chinese_mode_radiation(
    const ModalMode& mode, const ChineseRadiationConfig& config);

// Evaluates an already tracked candidate pair. Strike projections are mass-normalized
// modal participation magnitudes for the two branches (target, competitor).
ChineseBellObjectiveResult evaluate_chinese_bell_pair(
    const ModalMode& long_branch, const ModalMode& short_branch,
    double long_strike_on_long, double long_strike_on_short,
    double short_strike_on_long, double short_strike_on_short,
    const ChineseBellObjectiveConfig& config);

} // namespace BellOptimizer::Excitation
