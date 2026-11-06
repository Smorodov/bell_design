#include "listener_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Acoustics {
namespace {
ReceiverPoint makePoint(double x, double y, double z, double referenceDistance,
                        double soundSpeed) {
    ReceiverPoint p;
    p.x_m = x; p.y_m = y; p.z_m = z;
    p.distance_m = std::sqrt(x*x + y*y + z*z);
    if (!(p.distance_m > 0.0) || !std::isfinite(p.distance_m))
        throw std::invalid_argument("listener receiver is at the acoustic source origin");
    p.polar_rad = std::acos(std::clamp(z / p.distance_m, -1.0, 1.0));
    p.azimuth_rad = std::atan2(y, x);
    p.propagation_delay_s = (p.distance_m - referenceDistance) / soundSpeed;
    return p;
}
}

StereoListenerGeometry makeStereoListenerGeometry(double distance, double polar,
                                                    double azimuth, double spacing,
                                                    double yaw, double soundSpeed) {
    constexpr double halfPi = 1.57079632679489661923;
    if (!(std::isfinite(distance) && distance > 0.25) ||
        !(std::isfinite(polar) && polar >= 0.0 && polar <= 2.0 * halfPi) ||
        !std::isfinite(azimuth) ||
        !(std::isfinite(spacing) && spacing > 0.0 && spacing < 0.5) ||
        !(std::isfinite(yaw) && yaw >= -halfPi && yaw <= halfPi) ||
        !(std::isfinite(soundSpeed) && soundSpeed > 0.0))
        throw std::invalid_argument("invalid stereo listener geometry parameters");

    const double cx = distance * std::sin(polar) * std::cos(azimuth);
    const double cy = distance * std::sin(polar) * std::sin(azimuth);
    const double cz = distance * std::cos(polar);
    const double radialX = std::cos(azimuth), radialY = std::sin(azimuth);
    const double tangentX = -std::sin(azimuth), tangentY = std::cos(azimuth);
    const double baselineX = std::cos(yaw) * tangentX + std::sin(yaw) * radialX;
    const double baselineY = std::cos(yaw) * tangentY + std::sin(yaw) * radialY;
    const double halfSpacing = 0.5 * spacing;

    StereoListenerGeometry result;
    result.centre = makePoint(cx, cy, cz, distance, soundSpeed);
    result.left = makePoint(cx - halfSpacing * baselineX,
                            cy - halfSpacing * baselineY, cz, distance, soundSpeed);
    result.right = makePoint(cx + halfSpacing * baselineX,
                             cy + halfSpacing * baselineY, cz, distance, soundSpeed);
    const double dot = (result.left.x_m * result.right.x_m +
                        result.left.y_m * result.right.y_m +
                        result.left.z_m * result.right.z_m) /
                       (result.left.distance_m * result.right.distance_m);
    result.ear_angular_separation_rad = std::acos(std::clamp(dot, -1.0, 1.0));
    return result;
}

} // namespace BellOptimizer::Acoustics
