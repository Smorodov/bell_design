#include <cmath>
#include <iostream>
#include <vector>
#include "math_constants.hpp"

int main() {
    struct Case { double angular_velocity_impulse; double qd; double hz; };
    const std::vector<Case> cases = {
        {1.0, 1.0, 100.0},
        {0.75, -2.0, 460.3340261},
        {3.2, 0.125, 863.776078},
        {1.1, -0.4, 1754.943487}
    };

    double max_error = 0.0;
    for (const auto& c : cases) {
        const double omega = 2.0 * bell::kPi * c.hz;
        const double old_weight = std::abs(c.angular_velocity_impulse * c.qd);
        const double new_weight = old_weight * omega;
        const double expected = omega;
        const double ratio = old_weight > 0.0 ? new_weight / old_weight : 0.0;
        const double rel_error = std::abs(ratio / expected - 1.0);
        max_error = std::max(max_error, rel_error);
        std::cout << "hz=" << c.hz
                  << " old_weight=" << old_weight
                  << " pressure_weight=" << new_weight
                  << " new_over_old=" << ratio << '\n';
    }

    // Frequency scaling: keeping J_theta and Qd fixed, pressure weight is exactly
    // proportional to omega. This is the A/B contract used by the synth.
    const double w1 = 100.0;
    const double w2 = 400.0;
    const double ratio = (w2) / (w1);
    std::cout << "frequency_scaling_ratio=" << ratio << '\n';
    std::cout << "max_relative_contract_error=" << max_error << '\n';

    return (std::abs(ratio - 4.0) < 1e-12 && max_error < 1e-12) ? 0 : 1;
}
