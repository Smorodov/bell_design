#pragma once

#include <complex>
#include <vector>

namespace BellOptimizer::Audio {

struct ModalPressureTerm {
    double omega_rad_s = 0.0;
    double alpha_per_second = 0.0;
    double displacement_at_release = 0.0;
    double velocity_at_release = 0.0;
    std::complex<double> radiation_transfer{};
};

// Evaluate pressure at an absolute sample time. Propagation delay is represented
// solely as a time shift; callers must not also multiply the transfer by an
// additional exp(-i*omega*delay) phase factor.
double evaluateModalPressureAt(const std::vector<ModalPressureTerm>& modes,
                               double absoluteTimeSeconds,
                               double releaseTimeSeconds,
                               double propagationDelaySeconds);

} // namespace BellOptimizer::Audio
