#pragma once

#include "modal_types.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace BellOptimizer {

// Selection is deliberately performed on ModalData produced for the final,
// physically scaled geometry. Frequencies in this structure are absolute Hz.
struct PhysicalModalSelectionConfig {
    double minimum_frequency_hz = 0.0;
    double maximum_frequency_hz = 20000.0;
    double sample_rate_hz = 48000.0;
    double nyquist_margin = 0.98;
    // Relative to the strongest returned acoustic coefficient. Weak modes are
    // retained by default because they may contribute to attack and timbre.
    double weak_mode_relative_threshold = 1.0e-3;
    bool preserve_weak_modes = true;
    // Numerical integrity gates. Zero diagnostic values mean "not supplied".
    double maximum_eigen_residual = 1.0e-4;
    double maximum_mass_orthogonality_error = 1.0e-3;
};

struct PhysicalModalSelectionReport {
    std::vector<ModalMode> selected_modes;
    std::size_t input_count = 0;
    std::size_t outside_frequency_band = 0;
    std::size_t numerical_rejections = 0;
    std::size_t weak_modes = 0;
    std::size_t selected_count = 0;
    double effective_maximum_frequency_hz = 0.0;
    // Shape/family purity is not inferred from amplitude or eigensolver
    // residuals. It requires a calibrated reference-family/subspace metric.
    bool shape_purity_evaluated = false;
    std::string note;
};

PhysicalModalSelectionReport selectPhysicalSynthesisModes(
    const ModalData& physical_modal_data,
    const PhysicalModalSelectionConfig& config = {});

} // namespace BellOptimizer
