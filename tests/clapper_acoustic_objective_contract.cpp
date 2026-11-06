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
    ClapperModalState state;
    const double prime = targets.prime_freq_target;

    // Hum is a degenerate two-mode subspace. Its combined pressure norm is 1,
    // independent of the arbitrary orthonormal basis returned by the solver.
    append_mode(state, prime * targets.freq_ratios[0] * (1.0 - 1.0e-5), 1.0 / std::sqrt(2.0), 0.0);
    append_mode(state, prime * targets.freq_ratios[0] * (1.0 + 1.0e-5), 0.0, 1.0 / std::sqrt(2.0));
    append_mode(state, prime, targets.amp_ratios[0], 0.0);
    append_mode(state, prime * targets.freq_ratios[1], targets.amp_ratios[1], 0.0);
    append_mode(state, prime * targets.freq_ratios[2], targets.amp_ratios[2], 0.0);
    append_mode(state, prime * targets.freq_ratios[3], targets.amp_ratios[3], 0.0);

    ClapperAcousticObjectiveResult result;
    if (!evaluateClapperAcousticObjective(state, targets, result) ||
        !result.valid || result.total_error > 1.0e-12 || result.modal_amplitudes.size() != 5) {
        std::cerr << "target-matched grouped clapper objective should be zero\n";
        return 1;
    }

    // Rotate the degenerate pair's modal response into a different basis.
    // The sum-of-squares group norm, and therefore the score, must not change.
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

    state.excitation.modal_pressure_magnitude[3] = 0.0;
    state.excitation.modal_pressure_cosine[3] = 0.0;
    state.excitation.modal_pressure_sine[3] = 0.0;
    if (evaluateClapperAcousticObjective(state, targets, result)) {
        std::cerr << "zero target-family response must be rejected\n";
        return 4;
    }

    std::cout << "clapper_acoustic_objective_contract: PASS\n";
    return 0;
}
