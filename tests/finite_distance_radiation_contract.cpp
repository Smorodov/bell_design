#include "acoustics/finite_distance_radiation.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::Acoustics;
    constexpr double pi = 3.14159265358979323846;
    ModalMode mode;
    mode.frequency_hz = 620.0;
    mode.outer_surface_harmonics = {0, 1, -1};
    mode.outer_surface_fourier_radial = {0.21, 0.08, -0.03, 0.27, -0.04, 0.06, 0.31, 0.02, 0.04};
    mode.outer_surface_fourier_axial = {-0.11, 0.02, 0.05, -0.09, 0.03, -0.02, -0.07, 0.01, -0.04};
    const std::vector<Point2D> profile{{120.0, 0.0}, {180.0, 95.0}, {240.0, 210.0}};
    const double polar = 0.73, azimuth = -0.41;
    const double range = 2.0e6;
    ReceiverPoint receiver;
    receiver.x_m = range * std::sin(polar) * std::cos(azimuth);
    receiver.y_m = range * std::sin(polar) * std::sin(azimuth);
    receiver.z_m = range * std::cos(polar);
    receiver.distance_m = range;

    const auto coarseKernel = evaluateFiniteDistanceRadiation(mode, profile, receiver, 180);
    const auto nearKernel = evaluateFiniteDistanceRadiation(mode, profile, receiver, 720);
    const auto refinedKernel = evaluateFiniteDistanceRadiation(mode, profile, receiver, 1440);
    const double quadratureError = std::abs(nearKernel-refinedKernel) /
        std::max(1.0e-30, std::abs(refinedKernel));
    if (!(quadratureError < 2.0e-4)) {
        std::cerr << "azimuth quadrature convergence error=" << quadratureError << '\n';
        return 4;
    }

    // Independent far-field midpoint quadrature using the same surface patch
    // convention. At this range the finite-distance result must converge to it.
    const double k = 2.0 * pi * mode.frequency_hz / 343.0;
    const unsigned nPhi = 720;
    const double dPhi = 2.0 * pi / nPhi;
    std::complex<double> far{};
    auto fvalue = [&](const std::vector<double>& v, std::size_t j, double phi) {
        return v[j*3] + v[j*3+1]*std::cos(phi) + v[j*3+2]*std::sin(phi);
    };
    for (std::size_t j = 0; j + 1 < profile.size(); ++j) {
        const double dr = profile[j+1].r-profile[j].r;
        const double dz = profile[j+1].z-profile[j].z;
        const double dsMm = std::hypot(dr,dz);
        const double nr = dz/dsMm, nz = -dr/dsMm;
        const double r = 0.5*(profile[j].r+profile[j+1].r)*1e-3;
        const double z = 0.5*(profile[j].z+profile[j+1].z)*1e-3;
        const double ds = dsMm*1e-3;
        for (unsigned a=0; a<nPhi; ++a) {
            const double phi=dPhi*(a+0.5);
            const double ur=0.5*(fvalue(mode.outer_surface_fourier_radial,j,phi)+fvalue(mode.outer_surface_fourier_radial,j+1,phi));
            const double uz=0.5*(fvalue(mode.outer_surface_fourier_axial,j,phi)+fvalue(mode.outer_surface_fourier_axial,j+1,phi));
            const double ndisp=ur*nr+uz*nz;
            const double projection = r*std::sin(polar)*std::cos(phi-azimuth)+z*std::cos(polar);
            const double phase=-k*projection;
            far += ndisp*r*ds*dPhi*std::complex<double>(std::cos(phase),std::sin(phase));
        }
    }
    const double relative = std::abs(nearKernel-far)/std::max(1e-30,std::abs(far));
    if (!(relative < 2.0e-5)) {
        std::cerr << "far-field convergence error=" << relative << '\n';
        return 1;
    }

    // Moving the receiver while keeping the same ray direction changes phase
    // and amplitude continuously; invalid quadrature settings are rejected.
    receiver.x_m *= 0.5; receiver.y_m *= 0.5; receiver.z_m *= 0.5; receiver.distance_m *= 0.5;
    const auto closer = evaluateFiniteDistanceRadiation(mode, profile, receiver, 720);
    if (!std::isfinite(closer.real()) || !std::isfinite(closer.imag()) || std::abs(closer-nearKernel) < 1e-12 * std::max(1.0, std::abs(nearKernel))) return 2;
    bool rejected = false;
    try { (void)evaluateFiniteDistanceRadiation(mode, profile, receiver, 8); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) return 3;
    bool rejectedMeridional = false;
    try { (void)evaluateFiniteDistanceRadiation(mode, profile, receiver, 720, 343.0, 257); }
    catch (const std::invalid_argument&) { rejectedMeridional = true; }
    if (!rejectedMeridional) return 5;

    std::cout << "finite_distance_radiation_contract: PASS (far-field relative error=" << relative
              << ", quadrature relative error=" << quadratureError << ")\n";
    return 0;
}
