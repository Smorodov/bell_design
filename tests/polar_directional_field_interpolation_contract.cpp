#include "acoustics/directional_radiation.hpp"
#include "acoustics/listener_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::Acoustics;

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

ModalMode makeMode() {
    ModalMode mode;
    mode.frequency_hz = 913.0;
    mode.outer_surface_harmonics = {0, 1, -1, 2, -2, 3, -3};
    const std::size_t points = 6, basis = mode.outer_surface_harmonics.size();
    mode.outer_surface_fourier_radial.resize(points * basis);
    mode.outer_surface_fourier_axial.resize(points * basis);
    for (std::size_t p = 0; p < points; ++p) {
        for (std::size_t b = 0; b < basis; ++b) {
            const double x = 0.23 * static_cast<double>((p + 2) * (b + 1));
            mode.outer_surface_fourier_radial[p * basis + b] = 0.0007 * std::sin(x);
            mode.outer_surface_fourier_axial[p * basis + b] = 0.0004 * std::cos(0.71 * x);
        }
    }
    return mode;
}
}

int main() {
    try {
        const auto mode = makeMode();
        const std::vector<Point2D> profile{{35.0, 0.0}, {59.0, 88.0}, {91.0, 181.0},
                                           {128.0, 278.0}, {169.0, 379.0}, {205.0, 480.0}};
        const auto polarField = buildDirectionalRadiationPolarField(mode, profile, 2.0, 91);
        require(polarField.valid(), "polar directional field is invalid");

        double maxAbsolute = 0.0;
        unsigned comparisons = 0;
        double referencePeak = 0.0;
        for (double polar : {0.11, 0.37, 0.83, 1.17, 1.51, 2.03, 2.71, 3.02}) {
            const auto exactField = buildDirectionalRadiationField(mode, profile, polar, 2.0);
            for (double azimuth : {-2.2, -0.61, 0.0, 0.73, 2.41}) {
                const auto interpolated = polarField.evaluate(polar, azimuth);
                const auto exact = exactField.evaluate(azimuth);
                const double error = std::abs(interpolated - exact);
                referencePeak = std::max(referencePeak, std::abs(exact));
                maxAbsolute = std::max(maxAbsolute, error);
                ++comparisons;
            }
        }

        // The production renderer builds a local interval covering the complete
        // supported head-yaw range. Query multiple yaw values from that one field.
        double localMaxOverPeak = 0.0;
        double localReferencePeak = 0.0;
        for (double distance : {0.75, 1.5, 3.0, 6.0}) {
            const auto yawMinus = makeStereoListenerGeometry(distance, 0.83, -0.41, 0.18,
                -0.5 * 3.14159265358979323846);
            const auto yawPlus = makeStereoListenerGeometry(distance, 0.83, -0.41, 0.18,
                 0.5 * 3.14159265358979323846);
            const double pMin = std::min({yawMinus.left.polar_rad, yawMinus.right.polar_rad,
                                          yawPlus.left.polar_rad, yawPlus.right.polar_rad});
            const double pMax = std::max({yawMinus.left.polar_rad, yawMinus.right.polar_rad,
                                          yawPlus.left.polar_rad, yawPlus.right.polar_rad});
            const double span = pMax - pMin;
            const unsigned count = span < 1.0e-12 ? 1U :
                std::max(2U, static_cast<unsigned>(std::ceil(span / (2.0 * 3.14159265358979323846 / 180.0))) + 1U);
            const auto local = buildDirectionalRadiationPolarFieldRange(
                mode, profile, distance, pMin, pMax, count);
            for (double yaw : {-1.5707963267948966, -0.7853981633974483, 0.0,
                                0.7853981633974483, 1.5707963267948966}) {
                const auto ears = makeStereoListenerGeometry(distance, 0.83, -0.41, 0.18, yaw);
                for (const auto* ear : {&ears.left, &ears.right}) {
                    const auto exact = evaluateDirectionalRadiationDirect(
                        mode, profile, ear->polar_rad, ear->azimuth_rad, 4096);
                    const auto value = local.evaluate(ear->polar_rad, ear->azimuth_rad);
                    localReferencePeak = std::max(localReferencePeak, std::abs(exact));
                    localMaxOverPeak = std::max(localMaxOverPeak, std::abs(value - exact));
                }
            }
            const double midPolar = 0.5 * (pMin + pMax);
            const auto exactMid = evaluateDirectionalRadiationDirect(
                mode, profile, midPolar, -0.41, 4096);
            const auto interpolatedMid = local.evaluate(midPolar, -0.41);
            localReferencePeak = std::max(localReferencePeak, std::abs(exactMid));
            localMaxOverPeak = std::max(localMaxOverPeak, std::abs(interpolatedMid - exactMid));
        }
        std::cout << "local_yaw_envelope_error_over_peak="
                  << localMaxOverPeak / std::max(localReferencePeak, 1.0e-30) << '\n';
        require(localMaxOverPeak / std::max(localReferencePeak, 1.0e-30) < 0.005,
                "local polar interpolation exceeds 0.5% peak-normalized error across yaw envelope");

        const double rotation = 0.43;
        const double polar = 0.91, azimuth = 0.62;
        const auto rotated = polarField.evaluateRotated(polar, azimuth, rotation);
        const auto shifted = polarField.evaluate(polar, azimuth - rotation);
        require(std::abs(rotated - shifted) < 1.0e-13 * std::max(1.0, std::abs(shifted)),
                "polar interpolation and azimuthal frame rotation do not commute");

        std::cout << "polar_grid_samples=" << polarField.fields.size()
                  << " comparisons=" << comparisons
                  << " reference_peak=" << referencePeak
                  << " max_abs_over_peak=" << maxAbsolute / std::max(referencePeak, 1.0e-30)
                  << " local_yaw_envelope_error_over_peak=" << localMaxOverPeak / std::max(localReferencePeak, 1.0e-30) << '\n';
        require(maxAbsolute / std::max(referencePeak, 1.0e-30) < 0.005,
                "2-degree polar interpolation exceeds 0.5% peak-normalized complex error on fixture");
        std::cout << "polar_directional_field_interpolation_contract: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "polar_directional_field_interpolation_contract: FAIL: " << error.what() << '\n';
        return 1;
    }
}
