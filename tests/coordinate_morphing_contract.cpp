#include "math/coordinate_morphing.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using BellOptimizer::Math::AcousticWaveVector3D;
    using BellOptimizer::Math::CoordinateScale3D;

    const double pi = 3.14159265358979323846;
    const AcousticWaveVector3D k{2.0, 3.0, 6.0};

    // Identity must reproduce the original wave number exactly to floating-point
    // precision. This is the n=0/base-model gate for the morphing layer.
    const CoordinateScale3D identity{};
    const double k0 = BellOptimizer::Math::wave_number(k);
    const double k1 = BellOptimizer::Math::transformed_wave_number(k, identity);
    assert(std::abs(k0 - k1) < 1e-14);

    // Compressing y by 0.5 doubles only the y component of the wave vector.
    // This is the coordinate-pullback relation k'_y = k_y / s_y.
    const CoordinateScale3D compressed_y{1.0, 0.5, 1.0};
    const auto kp = BellOptimizer::Math::transform_wave_vector(k, compressed_y);
    assert(std::abs(kp.kx - 2.0) < 1e-14);
    assert(std::abs(kp.ky - 6.0) < 1e-14);
    assert(std::abs(kp.kz - 6.0) < 1e-14);

    const double expected = std::sqrt(2.0 * 2.0 + 6.0 * 6.0 + 6.0 * 6.0);
    assert(std::abs(BellOptimizer::Math::transformed_wave_number(k, compressed_y) - expected) < 1e-14);

    const double c = 343.0;
    const double f = BellOptimizer::Math::transformed_acoustic_frequency_hz(k, compressed_y, c);
    assert(std::abs(f - c * expected / (2.0 * pi)) < 1e-12);

    // A pure y-directed component demonstrates the expected inverse scaling:
    // lambda'_y = s_y lambda_y and therefore f'_y = f_y / s_y.
    const AcousticWaveVector3D ky{0.0, 4.0, 0.0};
    const double fy = BellOptimizer::Math::transformed_acoustic_frequency_hz(ky, compressed_y, c);
    const double fy_expected = c * 8.0 / (2.0 * pi);
    assert(std::abs(fy - fy_expected) < 1e-12);

    const auto point = compressed_y.point({10.0, 20.0, 30.0});
    assert(std::abs(point[0] - 10.0) < 1e-14);
    assert(std::abs(point[1] - 10.0) < 1e-14);
    assert(std::abs(point[2] - 30.0) < 1e-14);
    assert(std::abs(compressed_y.volume_scale() - 0.5) < 1e-14);

    std::cout << "coordinate morphing contract: PASS\n";
    return 0;
}
