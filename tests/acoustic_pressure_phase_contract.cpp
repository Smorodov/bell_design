#include <cmath>
#include <iostream>
#include "math_constants.hpp"

int main() {
    const double omega = 2.0 * bell::kPi * 863.776078;
    const double part_rz = -2.5;
    const double qd = 0.75;

    // Pressure coefficient for impulsive excitation and monopole radiation.
    const double c = -part_rz * qd * omega;
    const double flipped = -(-part_rz) * (-qd) * omega;

    // Eigenvectors have arbitrary sign. Both participation and Qd must flip,
    // leaving their physical product unchanged.
    if (!std::isfinite(c) || !std::isfinite(flipped)) return 1;
    if (std::abs(c - flipped) > 1e-12 * std::max(1.0, std::abs(c))) return 2;

    const double phase = (c < 0.0) ? bell::kPi : 0.0;
    if (!(phase == 0.0 || phase == bell::kPi)) return 3;

    std::cout << "pressure_coefficient=" << c << "\n";
    std::cout << "sign_flip_invariant_error=" << std::abs(c - flipped) << "\n";
    std::cout << "phase=" << phase << "\n";
    return 0;
}
