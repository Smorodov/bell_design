#include <cmath>
#include <iostream>
#include <limits>
#include "math_constants.hpp"

int main() {
    const double f1 = 460.0;
    const double f2 = 920.0;
    const double angular_velocity_impulse = 2.0;
    const double qd = 3.0;
    const double a1 = std::abs(angular_velocity_impulse * qd * (2.0 * bell::kPi * f1));
    const double a2 = std::abs(angular_velocity_impulse * qd * (2.0 * bell::kPi * f2));
    if (!std::isfinite(a1) || !std::isfinite(a2) || a1 <= 0.0) return 1;
    const double ratio = a2 / a1;
    if (std::abs(ratio - 2.0) > 1e-12) return 2;
    std::cout << "pressure_weight_ratio=" << ratio << "\n";
    return 0;
}
