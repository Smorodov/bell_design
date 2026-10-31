#include "chinese_bell_objective.hpp"

#include "math_constants.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {
namespace {
bool finite_positive(double x) { return std::isfinite(x) && x > 0.0; }
double ratio_or_zero(double target, double competitor) {
    return target > 0.0 ? competitor / target : std::numeric_limits<double>::infinity();
}
}

ChineseRadiationResult evaluate_chinese_mode_radiation(
    const ModalMode& mode, const ChineseRadiationConfig& c) {
    ChineseRadiationResult out;
    if (!finite_positive(mode.frequency_hz) || !finite_positive(c.air_density_kg_m3) ||
        !finite_positive(c.sound_speed_m_s) || !finite_positive(c.observation_distance_m) ||
        !std::isfinite(c.axial_direction_cosine) ||
        std::abs(c.axial_direction_cosine) > 1.0 ||
        !std::isfinite(mode.acoustic_volume_displacement) ||
        !std::isfinite(mode.acoustic_dipole_moment_z) ||
        !finite_positive(mode.decay_time)) return out;

    const double omega = 2.0 * bell::kPi * mode.frequency_hz;
    const double wave_number = omega / c.sound_speed_m_s;
    // The monopole and dipole are in quadrature in the far-field expansion.
    // Both terms have units m^3 per modal coordinate.
    out.monopole_m3 = mode.acoustic_volume_displacement;
    out.dipole_m3 = wave_number * c.axial_direction_cosine * mode.acoustic_dipole_moment_z;
    const double radiation_volume = std::hypot(out.monopole_m3, out.dipole_m3);
    out.pressure_pa_per_modal_coordinate =
        c.air_density_kg_m3 * omega * omega * radiation_volume /
        (4.0 * bell::kPi * c.observation_distance_m);
    out.t60_seconds = mode.decay_time;
    out.valid = std::isfinite(out.pressure_pa_per_modal_coordinate) &&
                out.pressure_pa_per_modal_coordinate >= 0.0;
    return out;
}

ChineseBellObjectiveResult evaluate_chinese_bell_pair(
    const ModalMode& a, const ModalMode& b,
    double a_on_a, double a_on_b, double b_on_a, double b_on_b,
    const ChineseBellObjectiveConfig& c) {
    ChineseBellObjectiveResult out;
    if (!finite_positive(a.frequency_hz) || !finite_positive(b.frequency_hz) ||
        !finite_positive(c.target_interval_cents) ||
        !finite_positive(c.interval_tolerance_cents) ||
        !finite_positive(c.minimum_t60_seconds) ||
        !finite_positive(c.maximum_competitor_ratio) ||
        !std::isfinite(a_on_a) || !std::isfinite(a_on_b) ||
        !std::isfinite(b_on_a) || !std::isfinite(b_on_b) ||
        a_on_a < 0.0 || a_on_b < 0.0 || b_on_a < 0.0 || b_on_b < 0.0) return out;

    out.interval_cents = 1200.0 * std::log2(b.frequency_hz / a.frequency_hz);
    out.interval_error_cents = std::abs(out.interval_cents - c.target_interval_cents);
    const auto ra = evaluate_chinese_mode_radiation(a, c.radiation);
    const auto rb = evaluate_chinese_mode_radiation(b, c.radiation);
    // An impulsive generalized load J gives modal displacement J/omega for
    // unit modal mass. Convert pressure-per-coordinate to pressure-per-unit
    // generalized impulse here; otherwise the strike amplitude is dimensionally
    // inconsistent and gains an erroneous extra factor of omega.
    out.long_strike_radiated_pressure =
        ra.pressure_pa_per_modal_coordinate * a_on_a / (2.0 * bell::kPi * a.frequency_hz);
    out.short_strike_radiated_pressure =
        rb.pressure_pa_per_modal_coordinate * b_on_b / (2.0 * bell::kPi * b.frequency_hz);
    out.long_branch_selectivity = ratio_or_zero(a_on_a, a_on_b);
    out.short_branch_selectivity = ratio_or_zero(b_on_b, b_on_a);
    out.weaker_target_t60_seconds = std::min(a.decay_time, b.decay_time);

    out.interval_gate = std::isfinite(out.interval_cents) &&
                        out.interval_error_cents <= c.interval_tolerance_cents;
    out.excitation_gate = a_on_a > 0.0 && b_on_b > 0.0 &&
        out.long_branch_selectivity <= c.maximum_competitor_ratio &&
        out.short_branch_selectivity <= c.maximum_competitor_ratio;
    out.radiation_gate = ra.valid && rb.valid &&
        out.long_strike_radiated_pressure > 0.0 &&
        out.short_strike_radiated_pressure > 0.0;
    out.decay_gate = finite_positive(out.weaker_target_t60_seconds) &&
                     out.weaker_target_t60_seconds >= c.minimum_t60_seconds;
    out.passed = out.interval_gate && out.excitation_gate &&
                 out.radiation_gate && out.decay_gate;
    return out;
}

} // namespace BellOptimizer::Excitation
