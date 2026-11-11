#include "core/physical_modal_selection.hpp"
#include <cstdlib>
#include <iostream>

using namespace BellOptimizer;

namespace {
void require(bool value, const char* message) {
    if (!value) { std::cerr << "FAIL: " << message << '\n'; std::exit(1); }
}
}

int main() {
    ModalData data;
    ModalMode low; low.frequency_hz = 440.0; low.acoustic_amplitude = 1.0;
    ModalMode weak; weak.frequency_hz = 880.0; weak.acoustic_amplitude = 1.0e-4;
    ModalMode high; high.frequency_hz = 23900.0; high.acoustic_amplitude = 0.5;
    ModalMode bad; bad.frequency_hz = 1200.0; bad.acoustic_amplitude = 0.4;
    bad.eigen_residual_relative = 0.1;
    data.modes = {high, weak, bad, low};

    PhysicalModalSelectionConfig config;
    config.sample_rate_hz = 48000.0;
    config.maximum_frequency_hz = 30000.0;
    config.weak_mode_relative_threshold = 1.0e-3;
    config.preserve_weak_modes = true;
    auto report = selectPhysicalSynthesisModes(data, config);
    require(report.effective_maximum_frequency_hz < 24000.0,
            "effective limit must respect Nyquist margin");
    require(report.selected_count == 2,
            "keep in-band modes and reject numerical outlier");
    require(report.weak_modes == 1,
            "weak mode must be identified independently of retention");
    require(report.selected_modes[0].frequency_hz == 440.0 &&
            report.selected_modes[1].frequency_hz == 880.0,
            "selection must preserve absolute-frequency order");
    require(!report.shape_purity_evaluated,
            "must not claim unimplemented modal shape purity");

    config.preserve_weak_modes = false;
    report = selectPhysicalSynthesisModes(data, config);
    require(report.selected_count == 1 && report.selected_modes[0].frequency_hz == 440.0,
            "weak-mode threshold must be optional and explicit");

    std::cout << "physical_modal_selection_contract: PASS\n";
    return 0;
}
