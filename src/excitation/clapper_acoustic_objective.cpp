#include "clapper_acoustic_objective.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace BellOptimizer::Excitation {
namespace {
bool finitePositive(double x) { return std::isfinite(x) && x > 0.0; }

// A perfectly circular bell has rotationally degenerate modal pairs. The
// eigensolver is free to return any orthonormal basis inside such a subspace;
// scoring the individual vectors would make the objective basis-dependent.
// Group only numerically coincident adjacent frequencies. The group pressure
// norm is invariant under an orthogonal rotation of an exactly degenerate pair.
constexpr double kRelativeDegeneracyTolerance = 5.0e-5;

struct ModalGroup {
    double frequency_hz = 0.0;
    double pressure_norm = 0.0;
    std::size_t count = 0;
};

bool build_modal_groups(const ClapperModalState& state, std::vector<ModalGroup>& groups) {
    groups.clear();
    const auto& modes = state.modal_data.modes;
    const auto& magnitudes = state.excitation.modal_pressure_magnitude;
    const auto& pressure_cosine = state.excitation.modal_pressure_cosine;
    const auto& pressure_sine = state.excitation.modal_pressure_sine;
    if (modes.empty() || magnitudes.size() != modes.size() ||
        pressure_cosine.size() != modes.size() || pressure_sine.size() != modes.size()) return false;

    double previous_frequency = 0.0;
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const double frequency = modes[i].frequency_hz;
        const double magnitude = magnitudes[i];
        if (!finitePositive(frequency) || frequency < previous_frequency ||
            !std::isfinite(magnitude) || magnitude < 0.0 ||
            !std::isfinite(pressure_cosine[i]) || !std::isfinite(pressure_sine[i])) return false;
        const double quadrature_magnitude = std::hypot(pressure_cosine[i], pressure_sine[i]);
        const double consistency_scale = std::max({1.0, magnitude, quadrature_magnitude});
        if (std::abs(quadrature_magnitude - magnitude) > 1.0e-9 * consistency_scale) return false;

        if (groups.empty()) {
            groups.push_back({frequency, magnitude * magnitude, 1});
            previous_frequency = frequency;
            continue;
        }

        ModalGroup& last = groups.back();
        const double relative_gap = std::abs(frequency - last.frequency_hz) /
            std::max(frequency, last.frequency_hz);
        if (relative_gap <= kRelativeDegeneracyTolerance) {
            const double n = static_cast<double>(last.count);
            last.frequency_hz = (last.frequency_hz * n + frequency) / (n + 1.0);
            last.pressure_norm += magnitude * magnitude;
            ++last.count;
        } else {
            groups.push_back({frequency, magnitude * magnitude, 1});
        }
        previous_frequency = frequency;
    }

    for (auto& group : groups) group.pressure_norm = std::sqrt(group.pressure_norm);
    return true;
}
} // namespace

bool evaluateClapperAcousticObjective(
    const ClapperModalState& state,
    const TargetData& targets,
    ClapperAcousticObjectiveResult& result) {
    result = {};
    std::vector<ModalGroup> groups;
    if (!build_modal_groups(state, groups) ||
        targets.frequency_reference_mode_index >= groups.size() ||
        !finitePositive(targets.prime_freq_target) ||
        (targets.frequency_targets.empty() && targets.amplitude_targets.empty())) return false;

    const double reference_frequency = groups[targets.frequency_reference_mode_index].frequency_hz;
    if (!finitePositive(reference_frequency)) return false;
    for (const auto& target : targets.frequency_targets) {
        if (target.weight == 0.0) continue;
        if (target.mode_index >= groups.size() || target.mode_index == targets.frequency_reference_mode_index ||
            !finitePositive(target.ratio) || !std::isfinite(target.weight) || target.weight < 0.0) return false;
        const double actual = groups[target.mode_index].frequency_hz / reference_frequency;
        const double error = std::log(actual / target.ratio);
        result.frequency_error += target.weight * error * error;
    }

    result.modal_amplitudes.resize(groups.size());
    double max_amplitude = 0.0;
    for (const auto& group : groups) max_amplitude = std::max(max_amplitude, group.pressure_norm);
    if (!finitePositive(max_amplitude)) return false;
    for (std::size_t i = 0; i < groups.size(); ++i)
        result.modal_amplitudes[i] = groups[i].pressure_norm / max_amplitude;

    for (const auto& target : targets.amplitude_targets) {
        if (target.weight == 0.0) continue;
        if (target.mode_index >= groups.size() || target.reference_mode_index >= groups.size() ||
            target.mode_index == target.reference_mode_index || !finitePositive(target.ratio) ||
            !std::isfinite(target.weight) || target.weight < 0.0) return false;
        const double numerator = groups[target.mode_index].pressure_norm;
        const double denominator = groups[target.reference_mode_index].pressure_norm;
        if (!finitePositive(numerator) || !finitePositive(denominator)) return false;
        const double actual = numerator / denominator;
        const double error = std::log(actual / target.ratio);
        result.amplitude_error += target.weight * error * error;
    }

    result.total_error = result.frequency_error + result.amplitude_error;
    result.valid = std::isfinite(result.total_error);
    return result.valid;
}

} // namespace BellOptimizer::Excitation
