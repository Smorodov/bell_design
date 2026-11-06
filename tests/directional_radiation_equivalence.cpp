#include "acoustics/directional_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        using BellOptimizer::ModalMode;
        using ::Point2D;
        using BellOptimizer::Acoustics::buildDirectionalRadiationField;
        using BellOptimizer::Acoustics::evaluateDirectionalRadiationDirect;

        ModalMode mode;
        mode.frequency_hz = 1873.25;
        mode.outer_surface_harmonics = {-3, -2, -1, 0, 1, 2, 3};
        const std::size_t points = 5;
        const std::size_t basis = mode.outer_surface_harmonics.size();
        mode.outer_surface_fourier_radial.resize(points * basis);
        mode.outer_surface_fourier_axial.resize(points * basis);
        for (std::size_t p = 0; p < points; ++p) {
            for (std::size_t b = 0; b < basis; ++b) {
                const double phase = 0.37 * static_cast<double>((p + 1) * (b + 2));
                mode.outer_surface_fourier_radial[p * basis + b] = 0.002 * std::sin(phase);
                mode.outer_surface_fourier_axial[p * basis + b] = 0.0015 * std::cos(phase * 0.83);
            }
        }
        const std::vector<Point2D> profile{{42.0, 0.0}, {68.0, 95.0}, {105.0, 210.0},
                                           {151.0, 330.0}, {196.0, 455.0}};
        double maxAbs = 0.0, maxRel = 0.0;
        for (double frequency : {320.0, 1873.25, 4200.0}) {
            mode.frequency_hz = frequency;
            for (double polar : {0.19, 0.71, 1.31}) {
                const auto field = buildDirectionalRadiationField(mode, profile, polar, 2.0);
                require(field.valid(), "harmonic directional field invalid");
                for (int k = 0; k < 17; ++k) {
                    const double azimuth = -1.1 + 0.173 * k;
                    const auto direct = evaluateDirectionalRadiationDirect(mode, profile, polar, azimuth, 8192);
                    const auto harmonic = field.evaluate(azimuth);
                    const double error = std::abs(direct - harmonic);
                    const double scale = std::max(1.0e-14, std::abs(direct));
                    maxAbs = std::max(maxAbs, error);
                    maxRel = std::max(maxRel, error / scale);
                }
            }
        }
        std::cout << "direct-vs-harmonic max_abs=" << maxAbs
                  << " max_rel=" << maxRel << '\n';
        require(maxRel < 2.0e-6, "direct and harmonic radiation disagree in complex amplitude/phase");
        const auto field = buildDirectionalRadiationField(mode, profile, 0.71, 2.0);

        // Rotation is an exact coordinate shift of the same retained field.
        const double baseAzimuth = 0.29, rotation = 0.41;
        const auto rotated = field.rotated(rotation);
        const auto expected = field.evaluate(baseAzimuth - rotation);
        const auto actual = rotated.evaluate(baseAzimuth);
        require(std::abs(expected - actual) < 1.0e-12 * std::max(1.0, std::abs(expected)),
                "rotated field does not preserve azimuthal evaluation");
        std::cout << "Directional radiation equivalence and rotation contracts passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Directional radiation equivalence failed: " << error.what() << '\n';
        return 1;
    }
}
