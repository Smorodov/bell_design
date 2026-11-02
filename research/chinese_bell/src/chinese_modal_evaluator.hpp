#pragma once

#include "fem/cylindrical_modal_evaluator.hpp"

namespace BellOptimizer::FEM {

struct ChineseBellStrikeSelectivityResult {
    std::vector<double> frequencies_hz;
    std::vector<std::size_t> mode_indices;
    std::vector<double> long_axis_excitation;
    std::vector<double> short_axis_excitation;
    double long_axis_cross_ratio = 0.0;
    double short_axis_cross_ratio = 0.0;
    bool valid = false;
};

struct ChineseBellStrikePositionResult {
    bool valid = false;
    std::size_t long_target_mode = 0;
    std::size_t short_target_mode = 0;
    double long_target_frequency_hz = 0.0;
    double short_target_frequency_hz = 0.0;
    double best_long_z = 0.0;
    double best_long_azimuth_rad = 0.0;
    double best_short_z = 0.0;
    double best_short_azimuth_rad = 0.0;
    double best_long_selectivity = 0.0;
    double best_short_selectivity = 0.0;
    double long_axis_selectivity = 0.0;
    double short_axis_selectivity = 0.0;
    double long_family_similarity = 0.0;
    double short_family_similarity = 0.0;
    bool family_tracking_valid = false;
    bool canonical_branch_separated = false;
    int canonical_branch_orientation = 0;
    double pair_frequency_split_cents = 0.0;
    double long_branch_contrast = 0.0;
    double short_branch_contrast = 0.0;
    bool two_tone_excitation_gate_passed = false;
};

bool evaluate_chinese_bell_strike_selectivity(
    const CylindricalBellModel& model,
    double center_z,
    double meridional_width,
    double angular_width_rad,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikeSelectivityResult& result);

bool optimize_chinese_bell_strike_pose(
    const CylindricalBellModel& model,
    double meridional_width,
    double angular_width_rad,
    unsigned z_samples,
    unsigned azimuth_samples,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikePositionResult& result,
    double canonical_strike_z_fraction = 0.70);

} // namespace BellOptimizer::FEM
