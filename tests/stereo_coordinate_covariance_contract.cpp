#include "acoustics/acoustic_directional_field.hpp"
#include "acoustics/listener_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
bool near(double a, double b, double tol = 2.0e-12) {
    return std::abs(a - b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}
bool nearComplex(const std::complex<double>& a, const std::complex<double>& b,
                 double tol = 2.0e-12) {
    return std::abs(a - b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}
struct Vec3 { double x, y, z; };
Vec3 rotateZ(double x, double y, double z, double angle) {
    return {x * std::cos(angle) - y * std::sin(angle),
            x * std::sin(angle) + y * std::cos(angle), z};
}
bool sameRotatedPoint(const BellOptimizer::Acoustics::ReceiverPoint& before,
                      const BellOptimizer::Acoustics::ReceiverPoint& after,
                      double angle) {
    const auto expected = rotateZ(before.x_m, before.y_m, before.z_m, angle);
    return near(expected.x, after.x_m) && near(expected.y, after.y_m) &&
           near(expected.z, after.z_m) && near(before.distance_m, after.distance_m) &&
           near(before.polar_rad, after.polar_rad) &&
           near(std::remainder(before.azimuth_rad + angle - after.azimuth_rad, 2.0*kPi), 0.0) &&
           near(before.propagation_delay_s, after.propagation_delay_s);
}
}

int main() {
    using namespace BellOptimizer::Acoustics;
    constexpr double radius = 2.7, polar = 0.63, azimuth = -0.42;
    constexpr double spacing = 0.18, yaw = 0.31, soundSpeed = 343.0;
    constexpr double sourceRotation = 1.137;

    const auto original = makeStereoListenerGeometry(radius, polar, azimuth,
                                                      spacing, yaw, soundSpeed);
    const auto rotated = makeStereoListenerGeometry(radius, polar,
        azimuth + sourceRotation, spacing, yaw, soundSpeed);
    if (!sameRotatedPoint(original.centre, rotated.centre, sourceRotation) ||
        !sameRotatedPoint(original.left, rotated.left, sourceRotation) ||
        !sameRotatedPoint(original.right, rotated.right, sourceRotation)) return 1;
    if (!near(original.ear_angular_separation_rad,
              rotated.ear_angular_separation_rad)) return 2;

    // A strike azimuth is a source-frame coordinate. Rotating the source,
    // strike and listener together must preserve their relative azimuth.
    const double strikeAzimuth = 0.74;
    const double listenerRelativeToStrike = azimuth - strikeAzimuth;
    const double rotatedStrikeAzimuth = strikeAzimuth + sourceRotation;
    const double rotatedRelative = (azimuth + sourceRotation) - rotatedStrikeAzimuth;
    if (!near(std::remainder(listenerRelativeToStrike - rotatedRelative, 2.0*kPi), 0.0)) return 3;

    // The complex receiver transfer is invariant when the field's coordinate
    // frame and the receiver azimuth are rotated by the same amount.
    AcousticDirectionalField field;
    field.observation_radius_m = radius;
    field.observation_polar_rad = polar;
    field.orders = {-4, -2, -1, 0, 1, 3, 4};
    field.coefficients = {{0.03, -0.01}, {-0.02, 0.04}, {0.01, 0.07},
                          {0.12, -0.03}, {-0.08, 0.02}, {0.04, 0.06},
                          {-0.01, -0.02}};
    if (!field.valid()) return 4;
    const auto originalLeft = field.evaluate(original.left.azimuth_rad);
    const auto originalRight = field.evaluate(original.right.azimuth_rad);
    const auto rotatedLeft = field.evaluateRotated(rotated.left.azimuth_rad, sourceRotation);
    const auto rotatedRight = field.evaluateRotated(rotated.right.azimuth_rad, sourceRotation);
    if (!nearComplex(originalLeft, rotatedLeft) ||
        !nearComplex(originalRight, rotatedRight)) return 5;

    // Geometry-induced interaural delays and complex phase factors are also
    // unchanged by a common z-axis rotation (no double-counted delay).
    if (!near(original.left.propagation_delay_s, rotated.left.propagation_delay_s) ||
        !near(original.right.propagation_delay_s, rotated.right.propagation_delay_s)) return 6;

    std::cout << "stereo_coordinate_covariance_contract: PASS (source/strike/listener frame, complex ear transfers, delays)\n";
    return 0;
}
