#pragma once

namespace BellOptimizer::Excitation {

enum class ModalDecayModel { FixedT60, StructuralLossFactor };

struct ModalDecayParameters {
    double alpha_per_second = 0.0;
    double t60_seconds = 0.0;
};

// FixedT60: alpha=ln(1000)/T60.
// StructuralLossFactor: alpha ~= eta*omega/2 = pi*eta*f,
// valid as a lightly damped, frequency-independent structural loss-factor
// approximation. T60 is the time for the free amplitude envelope to reach 0.001.
bool calculateModalDecay(double frequency_hz, ModalDecayModel model,
                          double fixed_t60_seconds, double structural_loss_factor,
                          ModalDecayParameters& result);

} // namespace BellOptimizer::Excitation
