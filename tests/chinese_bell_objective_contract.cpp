#include "excitation/chinese_bell_objective.hpp"
#include <cmath>
#include <iostream>

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::Excitation;
    ModalMode a, b;
    a.frequency_hz = 400.0; b.frequency_hz = 800.0;
    a.acoustic_volume_displacement = 1.0e-6;
    b.acoustic_volume_displacement = 0.8e-6;
    a.acoustic_dipole_moment_z = 2.0e-7;
    b.acoustic_dipole_moment_z = -1.0e-7;
    a.decay_time = 2.0; b.decay_time = 1.5;

    ChineseRadiationConfig rc;
    rc.axial_direction_cosine = 0.0;
    const auto side = evaluate_chinese_mode_radiation(a, rc);
    rc.axial_direction_cosine = 1.0;
    const auto axial = evaluate_chinese_mode_radiation(a, rc);
    if (!side.valid || !axial.valid || !(axial.pressure_pa_per_modal_coordinate > side.pressure_pa_per_modal_coordinate)) return 1;
    rc.observation_distance_m = 2.0;
    const auto far = evaluate_chinese_mode_radiation(a, rc);
    if (!far.valid || std::abs(far.pressure_pa_per_modal_coordinate / axial.pressure_pa_per_modal_coordinate - 0.5) > 1e-12) return 2;

    ChineseBellObjectiveConfig cfg;
    cfg.target_interval_cents = 1200.0;
    cfg.interval_tolerance_cents = 1.0;
    cfg.maximum_competitor_ratio = 0.05;
    cfg.minimum_t60_seconds = 1.0;
    cfg.radiation.axial_direction_cosine = 0.0;
    auto result = evaluate_chinese_bell_pair(a, b, 1.0, 0.01, 0.01, 1.0, cfg);
    if (!result.passed || !result.interval_gate || !result.excitation_gate ||
        !result.radiation_gate || !result.decay_gate) return 3;
    result = evaluate_chinese_bell_pair(a, b, 1.0, 0.2, 0.01, 1.0, cfg);
    if (result.passed || result.excitation_gate) return 4;
    b.decay_time = 0.2;
    result = evaluate_chinese_bell_pair(a, b, 1.0, 0.01, 0.01, 1.0, cfg);
    if (result.passed || result.decay_gate) return 5;
    b.decay_time = 1.5; b.frequency_hz = 760.0;
    result = evaluate_chinese_bell_pair(a, b, 1.0, 0.01, 0.01, 1.0, cfg);
    if (result.passed || result.interval_gate) return 6;
    std::cout << "chinese_bell_objective_contract: PASS\n"
              << "octave_interval_cents=" << 1200.0 * std::log2(2.0) << '\n'
              << "side_pressure=" << side.pressure_pa_per_modal_coordinate << '\n'
              << "axial_pressure=" << axial.pressure_pa_per_modal_coordinate << '\n';
    return 0;
}
