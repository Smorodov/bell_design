#include "modal_decay.hpp"
#include <cmath>

namespace BellOptimizer::Excitation {
namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kLn1000 = 6.9077552789821370521;
}

bool calculateModalDecay(double frequency_hz, ModalDecayModel model,
                         double fixed_t60_seconds, double structural_loss_factor,
                         ModalDecayParameters& result) {
    result = {};
    if (!std::isfinite(frequency_hz) || frequency_hz <= 0.0) return false;
    double alpha = 0.0;
    switch (model) {
    case ModalDecayModel::FixedT60:
        if (!std::isfinite(fixed_t60_seconds) || fixed_t60_seconds <= 0.0) return false;
        alpha = kLn1000 / fixed_t60_seconds;
        break;
    case ModalDecayModel::StructuralLossFactor:
        // Keep this model in the lightly damped regime; eta >= 0.1 is
        // outside the intended approximation and is rejected explicitly.
        if (!std::isfinite(structural_loss_factor) || structural_loss_factor <= 0.0 ||
            structural_loss_factor >= 0.1) return false;
        alpha = kPi * structural_loss_factor * frequency_hz;
        break;
    default:
        return false;
    }
    const double t60 = kLn1000 / alpha;
    if (!std::isfinite(alpha) || alpha <= 0.0 || !std::isfinite(t60) || t60 <= 0.0 ||
        alpha >= 2.0 * kPi * frequency_hz) return false;
    result.alpha_per_second = alpha;
    result.t60_seconds = t60;
    return true;
}

} // namespace BellOptimizer::Excitation
