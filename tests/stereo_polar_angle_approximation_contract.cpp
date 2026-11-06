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
constexpr double kPi = 3.1415926535897932384626433832795;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

ModalMode makeRadiatingMode() {
    ModalMode m;
    m.frequency_hz = 640.0;
    m.outer_surface_harmonics = {0, 1, -1, 2, -2};
    const std::size_t points = 3, basis = m.outer_surface_harmonics.size();
    m.outer_surface_fourier_radial.resize(points*basis);
    m.outer_surface_fourier_axial.resize(points*basis);
    const double radial[3][5] = {
        {0.8e-3, 0.35e-3, -0.12e-3, 0.08e-3, 0.04e-3},
        {1.1e-3, 0.48e-3, -0.16e-3, 0.11e-3, 0.05e-3},
        {0.6e-3, 0.22e-3, -0.09e-3, 0.05e-3, 0.03e-3}
    };
    const double axial[3][5] = {
        {0.12e-3, -0.04e-3, 0.03e-3, 0.02e-3, -0.01e-3},
        {0.18e-3, -0.06e-3, 0.04e-3, 0.03e-3, -0.02e-3},
        {0.09e-3, -0.03e-3, 0.02e-3, 0.01e-3, -0.01e-3}
    };
    for (std::size_t p=0; p<points; ++p)
        for (std::size_t b=0; b<basis; ++b) {
            m.outer_surface_fourier_radial[p*basis+b] = radial[p][b];
            m.outer_surface_fourier_axial[p*basis+b] = axial[p][b];
        }
    return m;
}

double absError(const std::complex<double>& a, const std::complex<double>& b) {
    return std::abs(a-b);
}
}

int main() {
    try {
        // Meridional profile is in millimetres for the radiation integral.
        const std::vector<Point2D> outer{{300.0, 0.0}, {270.0, 520.0}, {220.0, 1000.0}};
        const auto mode = makeRadiatingMode();
        constexpr double centrePolar = 0.83;
        constexpr double centreAzimuth = -0.41;
        constexpr double earSpacing = 0.18;
        constexpr double yaw = 0.29;
        const double distances[] = {0.75, 1.5, 3.0, 6.0};
        double globalReferencePeak = 0.0;
        struct Sample { double distance; double approximateError; double exactBuilderError; double phaseError; };
        std::vector<Sample> samples;
        for (double distance : distances) {
            const auto ears = makeStereoListenerGeometry(distance, centrePolar, centreAzimuth,
                                                         earSpacing, yaw);
            const auto centreField = buildDirectionalRadiationField(mode, outer, centrePolar, distance);
            double maxApproxError = 0.0, maxBuilderError = 0.0, maxPhaseError = 0.0;
            std::vector<std::pair<std::complex<double>, std::complex<double>>> receiverValues;
            for (const auto* ear : {&ears.left, &ears.right}) {
                const auto approximate = centreField.evaluate(ear->azimuth_rad);
                const auto direct = evaluateDirectionalRadiationDirect(mode, outer,
                    ear->polar_rad, ear->azimuth_rad, 8192);
                const auto exactPolarField = buildDirectionalRadiationField(mode, outer,
                    ear->polar_rad, distance);
                const auto exactBuilder = exactPolarField.evaluate(ear->azimuth_rad);
                maxApproxError = std::max(maxApproxError, absError(approximate, direct));
                maxBuilderError = std::max(maxBuilderError, absError(exactBuilder, direct));
                receiverValues.emplace_back(approximate, direct);
                globalReferencePeak = std::max(globalReferencePeak, std::abs(direct));
            }
            double localReferencePeak = 0.0;
            for (const auto& v : receiverValues)
                localReferencePeak = std::max(localReferencePeak, std::abs(v.second));
            for (const auto& v : receiverValues) {
                if (std::abs(v.second) > 0.1 * localReferencePeak) {
                    const double phase = std::abs(std::remainder(std::arg(v.first)-std::arg(v.second), 2.0*kPi));
                    maxPhaseError = std::max(maxPhaseError, phase);
                }
            }
            samples.push_back({distance, maxApproxError, maxBuilderError, maxPhaseError});
        }
        require(globalReferencePeak > 0.0, "radiation fixture produced zero reference field");
        std::cout << "distance_m,centre_polar_approx_error_over_peak,polar_specific_builder_error_over_peak,max_phase_error_rad\n";
        double worstApprox = 0.0, worstBuilder = 0.0;
        for (const auto& s : samples) {
            const double e = s.approximateError / globalReferencePeak;
            const double b = s.exactBuilderError / globalReferencePeak;
            worstApprox = std::max(worstApprox, e);
            worstBuilder = std::max(worstBuilder, b);
            std::cout << s.distance << ',' << e << ',' << b << ',' << s.phaseError << '\n';
        }
        require(worstBuilder < 2.0e-11, "polar-specific Fourier/Bessel field disagrees with direct quadrature");
        require(samples.size() == 4, "unexpected receiver-distance sample count");
        for (std::size_t i=1; i<samples.size(); ++i)
            require(samples[i].approximateError < samples[i-1].approximateError,
                    "centre-polar approximation error did not decrease with receiver distance");
        require(worstApprox > 0.05, "fixture no longer exposes a measurable centre-polar approximation error");
        std::cout << "stereo_polar_angle_approximation_contract: PASS (centre-polar approximation measured; polar-specific field agrees with direct quadrature)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "stereo_polar_angle_approximation_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
