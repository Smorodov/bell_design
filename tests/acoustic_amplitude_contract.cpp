#include "modal_types.hpp"
#include "math_constants.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using BellOptimizer::ModalMode;

static double coeff(const ModalMode& m) {
    return -m.part_rz * m.acoustic_volume_displacement *
           (2.0 * bell::kPi * m.frequency_hz);
}

int main() {
    ModalMode a, b;
    a.frequency_hz = 460.0;
    a.part_rz = 0.31;
    a.acoustic_volume_displacement = 0.004;
    b.frequency_hz = 920.0;
    b.part_rz = -0.17;
    b.acoustic_volume_displacement = 0.006;

    const double ca = coeff(a);
    const double cb = coeff(b);
    if (!(std::isfinite(ca) && std::isfinite(cb) && std::abs(ca) > 0.0)) return 1;

    const double ratio = std::abs(cb / ca);
    const double expected = std::abs(
        (b.part_rz * b.acoustic_volume_displacement * b.frequency_hz) /
        (a.part_rz * a.acoustic_volume_displacement * a.frequency_hz));
    if (std::abs(ratio - expected) > 1e-12) return 2;

    // Eigenvector sign reversal changes both participation and Qd sign, so
    // the signed pressure coefficient and its magnitude must remain invariant.
    const double ca_flipped = -(-a.part_rz) * (-a.acoustic_volume_displacement) *
                              (2.0 * bell::kPi * a.frequency_hz);
    if (std::abs(ca_flipped - ca) > 1e-12 * std::max(1.0, std::abs(ca))) return 3;

    // The acoustic coefficient is linear in frequency for this displacement-
    // based modal impulse convention.
    const double old_f = a.frequency_hz;
    a.frequency_hz *= 2.0;
    if (std::abs(std::abs(coeff(a) / ca) - 2.0) > 1e-12) return 4;
    a.frequency_hz = old_f;

    std::cout << "acoustic_amplitude_contract: PASS\n";
    return 0;
}
