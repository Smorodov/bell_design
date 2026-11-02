#include "acoustics/acoustic_directional_field.hpp"

#include <cmath>
#include <complex>
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
    field.observation_elevation_rad = 0.25;
    field.orders = {0, 1, 2, 4};
    field.coefficients = {{0.3, 0.0}, {0.2, -0.1}, {-0.1, 0.25}, {0.05, 0.07}};
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
    if (std::abs(rotated.coefficients[0] - field.coefficients[0]) > 1e-14) return 5;
    const auto expected_m2 = field.coefficients[2] *
        std::polar(1.0, -2.0 * theta0);
    if (std::abs(rotated.coefficients[2] - expected_m2) > 1e-14) return 6;

    // Rotation preserves each coefficient magnitude.
    for (std::size_t i = 0; i < field.size(); ++i) {
        if (std::abs(std::abs(rotated.coefficients[i]) -
                     std::abs(field.coefficients[i])) > 1e-14) return 7;
    }

    std::cout << "acoustic_directional_field_rotation_contract: PASS\n";
    return 0;
}
