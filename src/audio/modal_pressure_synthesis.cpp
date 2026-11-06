#include "modal_pressure_synthesis.hpp"

#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Audio {

double evaluateModalPressureAt(const std::vector<ModalPressureTerm>& modes,
                               double absoluteTimeSeconds,
                               double releaseTimeSeconds,
                               double propagationDelaySeconds) {
    if (!std::isfinite(absoluteTimeSeconds) || !std::isfinite(releaseTimeSeconds) ||
        !std::isfinite(propagationDelaySeconds))
        throw std::invalid_argument("invalid modal pressure time or propagation delay");
    const double localTime = absoluteTimeSeconds - releaseTimeSeconds - propagationDelaySeconds;
    if (localTime < 0.0) return 0.0;
    double pressure = 0.0;
    for (const auto& mode : modes) {
        if (!(std::isfinite(mode.omega_rad_s) && mode.omega_rad_s > 0.0) ||
            !(std::isfinite(mode.alpha_per_second) && mode.alpha_per_second >= 0.0) ||
            mode.omega_rad_s * mode.omega_rad_s <= mode.alpha_per_second * mode.alpha_per_second ||
            !std::isfinite(mode.displacement_at_release) || !std::isfinite(mode.velocity_at_release) ||
            !std::isfinite(mode.radiation_transfer.real()) || !std::isfinite(mode.radiation_transfer.imag()))
            throw std::invalid_argument("invalid damped modal pressure term");
        const double wd = std::sqrt(mode.omega_rad_s * mode.omega_rad_s -
                                    mode.alpha_per_second * mode.alpha_per_second);
        const double b = (mode.velocity_at_release + mode.alpha_per_second *
                          mode.displacement_at_release) / wd;
        const std::complex<double> exponent(-mode.alpha_per_second, wd);
        const std::complex<double> qddBase(mode.displacement_at_release, -b);
        const std::complex<double> qdd = qddBase * exponent * exponent * std::exp(exponent * localTime);
        pressure += std::real(mode.radiation_transfer * qdd);
    }
    if (!std::isfinite(pressure)) throw std::runtime_error("non-finite modal acoustic pressure");
    return pressure;
}

} // namespace BellOptimizer::Audio
