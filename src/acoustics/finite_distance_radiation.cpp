#include "finite_distance_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Acoustics {
namespace {
constexpr double kPi = 3.1415926535897932384626433832795;

double fourierValue(const std::vector<double>& values, std::size_t point,
                    std::size_t pointCount, const std::vector<int>& labels,
                    double angle) {
    if (values.size() != pointCount * labels.size())
        throw std::invalid_argument("modal surface Fourier array has inconsistent dimensions");
    double result = 0.0;
    for (std::size_t b = 0; b < labels.size(); ++b) {
        const int order = labels[b];
        const double coefficient = values[point * labels.size() + b];
        if (order == 0) result += coefficient;
        else if (order > 0) result += coefficient * std::cos(order * angle);
        else result += coefficient * std::sin(-order * angle);
    }
    return result;
}
}

std::complex<double> evaluateFiniteDistanceRadiation(
        const ModalMode& mode, const std::vector<Point2D>& profile,
        const ReceiverPoint& receiver, unsigned azimuthSamples,
        double soundSpeed, unsigned meridionalSubdivisions) {
    if (profile.size() < 2 || mode.outer_surface_harmonics.empty() ||
        !(std::isfinite(mode.frequency_hz) && mode.frequency_hz > 0.0) ||
        !(std::isfinite(receiver.distance_m) && receiver.distance_m > 0.0) ||
        !(std::isfinite(receiver.x_m) && std::isfinite(receiver.y_m) && std::isfinite(receiver.z_m)) ||
        !(azimuthSamples >= 24 && azimuthSamples <= 16384) ||
        !(meridionalSubdivisions >= 1 && meridionalSubdivisions <= 256) ||
        !(std::isfinite(soundSpeed) && soundSpeed > 0.0))
        throw std::invalid_argument("invalid finite-distance acoustic evaluation input");
    if (mode.outer_surface_fourier_radial.size() != profile.size() * mode.outer_surface_harmonics.size() ||
        mode.outer_surface_fourier_axial.size() != profile.size() * mode.outer_surface_harmonics.size())
        throw std::invalid_argument("incomplete modal outer-surface Fourier field");

    const double k = 2.0 * kPi * mode.frequency_hz / soundSpeed;
    const double dPhi = 2.0 * kPi / static_cast<double>(azimuthSamples);
    const double receiverRadius = receiver.distance_m;
    std::complex<double> integral{};

    for (std::size_t j = 0; j + 1 < profile.size(); ++j) {
        const double drMm = profile[j + 1].r - profile[j].r;
        const double dzMm = profile[j + 1].z - profile[j].z;
        const double segmentLengthMm = std::hypot(drMm, dzMm);
        if (!(segmentLengthMm > 0.0)) continue;
        const double nr = dzMm / segmentLengthMm;
        const double nz = -drMm / segmentLengthMm;
        const double subLength = segmentLengthMm * 1.0e-3 /
                                 static_cast<double>(meridionalSubdivisions);

        // Composite midpoint integration along each existing meridional
        // profile segment. Geometry and modal displacement are linearly
        // interpolated; the segment normal is constant on this piecewise-
        // linear profile. Refinement therefore controls integration error,
        // not the underlying geometric/FEM profile resolution.
        for (unsigned m = 0; m < meridionalSubdivisions; ++m) {
            const double t = (static_cast<double>(m) + 0.5) /
                             static_cast<double>(meridionalSubdivisions);
            const double r = ((1.0 - t) * profile[j].r + t * profile[j + 1].r) * 1.0e-3;
            const double z = ((1.0 - t) * profile[j].z + t * profile[j + 1].z) * 1.0e-3;
            if (!(r >= 0.0) || !std::isfinite(r) || !std::isfinite(z))
                throw std::invalid_argument("invalid meridional profile coordinate");

            for (unsigned a = 0; a < azimuthSamples; ++a) {
                const double phi = dPhi * (static_cast<double>(a) + 0.5);
                const double ur0 = fourierValue(mode.outer_surface_fourier_radial, j,
                    profile.size(), mode.outer_surface_harmonics, phi);
                const double ur1 = fourierValue(mode.outer_surface_fourier_radial, j + 1,
                    profile.size(), mode.outer_surface_harmonics, phi);
                const double uz0 = fourierValue(mode.outer_surface_fourier_axial, j,
                    profile.size(), mode.outer_surface_harmonics, phi);
                const double uz1 = fourierValue(mode.outer_surface_fourier_axial, j + 1,
                    profile.size(), mode.outer_surface_harmonics, phi);
                const double ur = (1.0 - t) * ur0 + t * ur1;
                const double uz = (1.0 - t) * uz0 + t * uz1;
                const double normalDisplacement = ur * nr + uz * nz;
                const double x = r * std::cos(phi);
                const double y = r * std::sin(phi);
                const double dx = receiver.x_m - x;
                const double dy = receiver.y_m - y;
                const double dzReceiver = receiver.z_m - z;
                const double distance = std::sqrt(dx*dx + dy*dy + dzReceiver*dzReceiver);
                if (!(distance > 1.0e-9) || !std::isfinite(distance))
                    throw std::invalid_argument("receiver intersects the acoustic surface quadrature point");

                // exp(+ik(R-R0))*R0/R preserves the existing far-field phase
                // convention; R0 is the listener-center reference distance.
                const double amplitude = receiverRadius / distance;
                const double phase = k * (distance - receiverRadius);
                integral += normalDisplacement * r * subLength * dPhi * amplitude *
                            std::complex<double>(std::cos(phase), std::sin(phase));
            }
        }
    }
    if (!std::isfinite(integral.real()) || !std::isfinite(integral.imag()))
        throw std::runtime_error("non-finite finite-distance acoustic integral");
    return integral;
}

} // namespace BellOptimizer::Acoustics
