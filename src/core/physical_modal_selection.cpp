#include "physical_modal_selection.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer {

PhysicalModalSelectionReport selectPhysicalSynthesisModes(
    const ModalData& physical_modal_data,
    const PhysicalModalSelectionConfig& config) {
    if (!std::isfinite(config.minimum_frequency_hz) ||
        !std::isfinite(config.maximum_frequency_hz) ||
        !std::isfinite(config.sample_rate_hz) ||
        !std::isfinite(config.nyquist_margin) ||
        !std::isfinite(config.weak_mode_relative_threshold) ||
        config.minimum_frequency_hz < 0.0 ||
        config.maximum_frequency_hz <= config.minimum_frequency_hz ||
        config.sample_rate_hz <= 0.0 || config.nyquist_margin <= 0.0 ||
        config.nyquist_margin > 1.0 || config.weak_mode_relative_threshold < 0.0 ||
        config.maximum_eigen_residual < 0.0 ||
        config.maximum_mass_orthogonality_error < 0.0) {
        throw std::invalid_argument("invalid physical modal selection configuration");
    }

    PhysicalModalSelectionReport report;
    report.input_count = physical_modal_data.modes.size();
    report.effective_maximum_frequency_hz = std::min(
        config.maximum_frequency_hz,
        0.5 * config.sample_rate_hz * config.nyquist_margin);

    double maximum_acoustic = 0.0;
    for (const auto& mode : physical_modal_data.modes) {
        const double value = mode.acoustic_amplitude > 0.0
            ? mode.acoustic_amplitude : mode.amplitude;
        if (std::isfinite(value)) maximum_acoustic = std::max(maximum_acoustic, value);
    }

    for (const auto& mode : physical_modal_data.modes) {
        if (!std::isfinite(mode.frequency_hz) ||
            mode.frequency_hz < config.minimum_frequency_hz ||
            mode.frequency_hz > report.effective_maximum_frequency_hz) {
            ++report.outside_frequency_band;
            continue;
        }
        // A zero diagnostic means the solver/exporter did not provide it; do
        // not reject legacy modes on missing diagnostics alone.
        if (!std::isfinite(mode.eigen_residual_relative) ||
            !std::isfinite(mode.mass_orthogonality_error) ||
            (mode.eigen_residual_relative > 0.0 &&
             mode.eigen_residual_relative > config.maximum_eigen_residual) ||
            (mode.mass_orthogonality_error > 0.0 &&
             mode.mass_orthogonality_error > config.maximum_mass_orthogonality_error)) {
            ++report.numerical_rejections;
            continue;
        }
        const double acoustic = mode.acoustic_amplitude > 0.0
            ? mode.acoustic_amplitude : mode.amplitude;
        const double relative = maximum_acoustic > 0.0 ? acoustic / maximum_acoustic : 0.0;
        const bool weak = relative < config.weak_mode_relative_threshold;
        if (weak) ++report.weak_modes;
        if (weak && !config.preserve_weak_modes) continue;
        report.selected_modes.push_back(mode);
    }

    std::stable_sort(report.selected_modes.begin(), report.selected_modes.end(),
        [](const ModalMode& a, const ModalMode& b) {
            return a.frequency_hz < b.frequency_hz;
        });
    report.selected_count = report.selected_modes.size();
    report.shape_purity_evaluated = false;
    report.note = "Absolute-frequency and numerical-integrity selection applied after physical scaling; "
                  "modal shape/family purity requires a separately calibrated subspace metric.";
    return report;
}

} // namespace BellOptimizer
