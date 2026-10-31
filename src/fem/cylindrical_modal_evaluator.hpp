#pragma once

#include "cylindrical_bell_model.hpp"
#include "core/modal_types.hpp"
#include "excitation/chinese_bell_strike.hpp"

namespace BellOptimizer::FEM {

struct ChineseBellStrikeSelectivityResult {
    std::vector<double> frequencies_hz;
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
    // Mechanical two-tone excitation gate only. This does not imply that the
    // full Chinese acoustic objective is ready for production optimization.
    bool two_tone_excitation_gate_passed = false;
};

struct CylindricalModalEvaluationConfig {
    std::vector<int> harmonics{0, 1, 2, 3};
    unsigned angular_samples = 64;
    int requested_modes = 10;
    int max_iterations = 1200;
    double tolerance = 1.0e-7;
    // The optimizer canonicalizes geometry to H=1000 mm before direct FEM.
    // This gate removes residual rigid-body numerical modes below the production
    // minimum frequency without hiding physical modes in the canonical search space.
    double minimum_frequency_hz = 1.0;
    bool use_rcm_ordering = true;
};

// Direct production-side modal evaluation for the vector cylindrical FEM path.
// It deliberately contains no coordinate morphing or frequency post-correction.
// The returned amplitude is a mass-projected meridional influence diagnostic;
// acoustic radiation is kept outside this bridge until the Chinese-bell objective
// defines its strike and radiation contract explicitly.
bool evaluate_cylindrical_bell_modes(
    const CylindricalBellModel& model,
    const ModalInfluence& influence,
    const CylindricalModalEvaluationConfig& config,
    ModalData& result);

// Evaluate two canonical Chinese-bell strike positions: azimuth 0 (long axis)
// and pi/2 (short axis). The result is a modal-selectivity diagnostic; it does
// not yet alter the production optimizer objective.
bool evaluate_chinese_bell_strike_selectivity(
    const CylindricalBellModel& model,
    double center_z,
    double meridional_width,
    double angular_width_rad,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikeSelectivityResult& result);

// Search the outer wall for the best strike positions for the two modal
// responses identified by the canonical long/short-axis strikes. The FEM solve
// is performed once; only the modal load projection is swept over (z, azimuth).
bool optimize_chinese_bell_strike_pose(
    const CylindricalBellModel& model,
    double meridional_width,
    double angular_width_rad,
    unsigned z_samples,
    unsigned azimuth_samples,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikePositionResult& result);

} // namespace BellOptimizer::FEM
