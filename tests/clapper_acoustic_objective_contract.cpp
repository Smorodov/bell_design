#include "excitation/clapper_acoustic_objective.hpp"

#include <cmath>
#include <iostream>

using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;

static void append_mode(ClapperModalState& state, double frequency, double cosine, double sine) {
    ModalMode mode;
    mode.frequency_hz = frequency;
    state.modal_data.modes.push_back(mode);
    state.excitation.modal_pressure_cosine.push_back(cosine);
    state.excitation.modal_pressure_sine.push_back(sine);
    state.excitation.modal_pressure_magnitude.push_back(std::hypot(cosine, sine));
}

int main() {
    TargetData targets;
    targets.frequency_reference_mode_index = 1;
    targets.prime_freq_target = 440.0;
    targets.frequency_targets = {
        {0, 0.5, 1.0, "hum"}, {2, 1.25, 1.0, "third"},
        {3, 1.5, 1.0, "fifth"}, {4, 2.0, 1.0, "nominal"}
    };
    targets.amplitude_targets = {
        {1, 0, 0.8, 1.0, "prime"}, {2, 0, 0.6, 1.0, "third"},
        {3, 0, 0.4, 1.0, "fifth"}, {4, 0, 0.2, 1.0, "nominal"}
    };

    ClapperModalState state;
    const double prime = targets.prime_freq_target;
    // The hum is a degenerate two-mode subspace. The evaluator combines its
    // response by RSS before applying profile selectors.
    append_mode(state, prime * 0.5 * (1.0 - 1.0e-5), 1.0 / std::sqrt(2.0), 0.0);
    append_mode(state, prime * 0.5 * (1.0 + 1.0e-5), 0.0, 1.0 / std::sqrt(2.0));
    append_mode(state, prime, 0.8, 0.0);
    append_mode(state, prime * 1.25, 0.6, 0.0);
    append_mode(state, prime * 1.5, 0.4, 0.0);
    append_mode(state, prime * 2.0, 0.2, 0.0);

    ClapperAcousticObjectiveResult result;
    if (!evaluateClapperAcousticObjective(state, targets, result) ||
        !result.valid || result.total_error > 1.0e-12 || result.modal_amplitudes.size() != 5) {
        std::cerr << "target-matched grouped clapper objective should be zero\n";
        return 1;
    }

    // Rotate response inside the degenerate pair; group RSS and score stay fixed.
    state.excitation.modal_pressure_cosine[0] = 1.0;
    state.excitation.modal_pressure_sine[0] = 0.0;
    state.excitation.modal_pressure_magnitude[0] = 1.0;
    state.excitation.modal_pressure_cosine[1] = 0.0;
    state.excitation.modal_pressure_sine[1] = 0.0;
    state.excitation.modal_pressure_magnitude[1] = 0.0;
    if (!evaluateClapperAcousticObjective(state, targets, result) || result.total_error > 1.0e-12) {
        std::cerr << "degenerate-subspace objective must be basis invariant\n";
        return 2;
    }

    for (double& x : state.excitation.modal_pressure_magnitude) x *= 17.0;
    for (double& x : state.excitation.modal_pressure_cosine) x *= 17.0;
    for (double& x : state.excitation.modal_pressure_sine) x *= 17.0;
    if (!evaluateClapperAcousticObjective(state, targets, result) || result.total_error > 1.0e-12) {
        std::cerr << "objective must be invariant to common acoustic gain\n";
        return 3;
    }

    // A profile with a different target count and a non-default frequency
    // reference is evaluated by the same generic objective.
    TargetData compact;
    compact.frequency_reference_mode_index = 2;
    compact.prime_freq_target = 550.0;
    compact.frequency_targets = {{0, 0.4, 1.0, "low"}, {4, 1.6, 0.5, "high"}};
    if (!evaluateClapperAcousticObjective(state, compact, result) || result.total_error > 1.0e-12) {
        std::cerr << "variable target count/arbitrary reference must be supported\n";
        return 4;
    }

    targets.frequency_targets.push_back({99, 3.0, 1.0, "out_of_range"});
    if (evaluateClapperAcousticObjective(state, targets, result)) {
        std::cerr << "out-of-range generic selector must be rejected\n";
        return 5;
    }
    std::cout << "clapper_acoustic_objective_contract: PASS\n";
    return 0;
}
