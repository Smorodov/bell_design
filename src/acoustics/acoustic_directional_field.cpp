#include "acoustic_directional_field.hpp"

#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Acoustics {
namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kTwoPi = 2.0 * kPi;
}

bool AcousticDirectionalField::valid() const {
    if (!(std::isfinite(observation_radius_m) && observation_radius_m > 0.0)) return false;
    if (!std::isfinite(observation_elevation_rad)) return false;
    if (orders.size() != coefficients.size() || orders.empty()) return false;
    for (std::size_t i = 0; i < orders.size(); ++i) {
        if (orders[i] < 0 || !std::isfinite(coefficients[i].real()) ||
            !std::isfinite(coefficients[i].imag())) return false;
        if (i > 0 && orders[i] <= orders[i - 1]) return false;
    }
    return true;
}

AcousticDirectionalField AcousticDirectionalField::rotated(double theta_rad) const {
    if (!valid() || !std::isfinite(theta_rad))
        throw std::invalid_argument("invalid acoustic directional field or rotation");

    AcousticDirectionalField result = *this;
    const double wrapped = std::remainder(theta_rad, kTwoPi);
    for (std::size_t i = 0; i < result.coefficients.size(); ++i) {
        const double phase = -static_cast<double>(result.orders[i]) * wrapped;
        result.coefficients[i] *= std::polar(1.0, phase);
    }
    return result;
}

std::complex<double> AcousticDirectionalField::evaluate(double theta_rad) const {
    if (!valid() || !std::isfinite(theta_rad))
        throw std::invalid_argument("invalid acoustic directional field or azimuth");

    const double wrapped = std::remainder(theta_rad, kTwoPi);
    std::complex<double> value{0.0, 0.0};
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double phase = static_cast<double>(orders[i]) * wrapped;
        value += coefficients[i] * std::polar(1.0, phase);
    }
    return value;
}

} // namespace BellOptimizer::Acoustics
