#include "acoustics/acoustic_directional_field.hpp"

#include <cmath>
#include <complex>
#include <algorithm>
#include <iostream>

namespace {
constexpr double kPi = 3.1415926535897932384626433832795;

double error(const std::complex<double>& a, const std::complex<double>& b) {
    return std::abs(a - b) / std::max(1.0, std::abs(b));
}
}

int main() {
    using BellOptimizer::Acoustics::AcousticDirectionalField;
    AcousticDirectionalField field;
    field.observation_radius_m = 10.0;
    field.observation_polar_rad = 0.25;
    field.orders = {-3, -1, 0, 1, 2, 4};
    field.coefficients = {{0.02, -0.04}, {0.1, 0.07}, {0.3, 0.0}, {0.2, -0.1}, {-0.1, 0.25}, {0.05, 0.07}};
    if (!field.valid()) return 1;

    // Zero and full-turn rotations must leave the physical field unchanged.
    if (error(field.rotated(0.0).evaluate(0.37), field.evaluate(0.37)) > 1e-14) return 2;
    if (error(field.rotated(2.0 * kPi).evaluate(0.37), field.evaluate(0.37)) > 1e-14) return 3;

    // Rotation must be equivalent to evaluating the base field at theta-theta0.
    const double theta0 = 0.61;
    const double theta = -0.27;
    const auto rotated_value = field.rotated(theta0).evaluate(theta);
    const auto expected = field.evaluate(theta - theta0);
    if (error(rotated_value, expected) > 1e-14) return 4;

    // m=0 is invariant; nonzero orders must acquire their corresponding phase.
    const auto rotated = field.rotated(theta0);
    const auto zeroIndex = static_cast<std::size_t>(std::find(field.orders.begin(), field.orders.end(), 0) - field.orders.begin());
    const auto m2Index = static_cast<std::size_t>(std::find(field.orders.begin(), field.orders.end(), 2) - field.orders.begin());
    if (std::abs(rotated.coefficients[zeroIndex] - field.coefficients[zeroIndex]) > 1e-14) return 5;
    const auto expected_m2 = field.coefficients[m2Index] *
        std::polar(1.0, -2.0 * theta0);
    if (std::abs(rotated.coefficients[m2Index] - expected_m2) > 1e-14) return 6;

    // Rotation preserves each coefficient magnitude.
    for (std::size_t i = 0; i < field.size(); ++i) {
        if (std::abs(std::abs(rotated.coefficients[i]) -
                     std::abs(field.coefficients[i])) > 1e-14) return 7;
    }

    // Joint source/listener frame rotation is covariant: rotating both
    // coordinate frames leaves the complex transfer unchanged.
    const double frameRotation = -0.83;
    const double listenerAzimuth = 1.17;
    const auto covariant = field.rotated(frameRotation).evaluate(
        listenerAzimuth + frameRotation);
    if (error(covariant, field.evaluate(listenerAzimuth)) > 1e-14) return 8;

    // Repeated receiver queries use the immutable base coefficients. The
    // no-copy API must match a materialized rotation for arbitrary angles.
    const auto originalOrders = field.orders;
    const auto originalCoefficients = field.coefficients;
    for (double rotation : {-2.4, -0.3, 0.0, 0.91, 2.8}) {
        for (double receiver : {-2.1, -0.2, 0.6, 2.3}) {
            const auto noCopy = field.evaluateRotated(receiver, rotation);
            const auto materialized = field.rotated(rotation).evaluate(receiver);
            if (error(noCopy, materialized) > 1e-14) return 9;
        }
    }
    if (field.orders != originalOrders ||
        field.coefficients != originalCoefficients) return 10;

    // Rotations compose additively without changing the underlying field.
    const double a = 0.37, b = -1.12, thetaCheck = 0.49;
    const auto composed = field.evaluateRotated(thetaCheck, a + b);
    const auto sequential = field.rotated(a).evaluateRotated(thetaCheck, b);
    if (error(composed, sequential) > 1e-14) return 11;

    std::cout << "acoustic_directional_field_rotation_contract: PASS (covariance, composition, immutable reuse)\n";
    return 0;
}
