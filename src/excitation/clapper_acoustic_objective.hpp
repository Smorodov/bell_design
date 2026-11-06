#pragma once

#include "core/bell_types.h"
#include "excitation/clapper_modal_state_provider.hpp"

#include <vector>

namespace BellOptimizer::Excitation {

struct ClapperAcousticObjectiveResult {
    bool valid = false;
    double frequency_error = 0.0;
    double amplitude_error = 0.0;
    double total_error = 1.0e9;
    // Normalized root-sum-square amplitudes, one value per modal group.
    std::vector<double> modal_amplitudes;
};

// Evaluates modal-family frequency ratios and clapper-driven acoustic pressure.
// Numerically degenerate adjacent modes are grouped into one subspace: frequencies
// are averaged and pressure magnitudes are combined by root-sum-square, avoiding
// dependence on the arbitrary eigenvector basis inside a circular-bell pair. A
// common gain cancels in amplitude ratios; legacy ModalMode::amplitude is unused.
bool evaluateClapperAcousticObjective(
    const ClapperModalState& state,
    const TargetData& targets,
    ClapperAcousticObjectiveResult& result);

} // namespace BellOptimizer::Excitation
