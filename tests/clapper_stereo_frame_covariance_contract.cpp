#include "excitation/clapper_dynamics.hpp"
#include "acoustics/directional_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;
using namespace BellOptimizer::Acoustics;

namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
bool near(double a, double b, double tol = 3.0e-8) {
    return std::abs(a-b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}
bool nearComplex(std::complex<double> a, std::complex<double> b, double tol = 2.0e-8) {
    return std::abs(a-b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}

std::vector<ModalMode> makeBaseModes() {
    std::vector<ModalMode> modes(3);
    const double innerR[2] = {0.0010, 0.0012};
    const double innerZ[2] = {0.0001, 0.0002};
    const double outerR[2] = {0.0008, 0.0010};
    const double outerZ[2] = {0.0002, 0.0003};
    for (int k=0; k<2; ++k) {
        auto& m = modes[static_cast<std::size_t>(k)];
        m.frequency_hz = 180.0;
        m.modal_mass = 1.0;
        m.inner_surface_harmonics = {k == 0 ? 1 : -1};
        m.inner_surface_fourier_radial.assign(innerR, innerR+2);
        m.inner_surface_fourier_axial.assign(innerZ, innerZ+2);
        m.outer_surface_harmonics = {k == 0 ? 1 : -1};
        m.outer_surface_fourier_radial.assign(outerR, outerR+2);
        m.outer_surface_fourier_axial.assign(outerZ, outerZ+2);
    }
    auto& axis = modes[2];
    axis.frequency_hz = 410.0;
    axis.modal_mass = 1.0;
    axis.inner_surface_harmonics = {0};
    axis.inner_surface_fourier_radial = {0.0007, 0.0009};
    axis.inner_surface_fourier_axial = {0.0001, 0.0002};
    axis.outer_surface_harmonics = {0};
    axis.outer_surface_fourier_radial = {0.0005, 0.0007};
    axis.outer_surface_fourier_axial = {0.0002, 0.0003};
    return modes;
}

std::vector<ModalMode> rotateDegeneratePair(const std::vector<ModalMode>& base, double a) {
    auto out = base;
    const double c = std::cos(a), s = std::sin(a);
    for (std::size_t k=0; k<2; ++k) {
        const auto& sourceCos = base[0];
        const auto& sourceSin = base[1];
        auto& target = out[k];
        target.inner_surface_harmonics = {1, -1};
        target.outer_surface_harmonics = {1, -1};
        target.inner_surface_fourier_radial.assign(4, 0.0);
        target.inner_surface_fourier_axial.assign(4, 0.0);
        target.outer_surface_fourier_radial.assign(4, 0.0);
        target.outer_surface_fourier_axial.assign(4, 0.0);
        const double cc = (k == 0) ? c : -s;
        const double ss = (k == 0) ? s : c;
        for (std::size_t p=0; p<2; ++p) {
            target.inner_surface_fourier_radial[p*2] = cc * sourceCos.inner_surface_fourier_radial[p];
            target.inner_surface_fourier_radial[p*2+1] = ss * sourceSin.inner_surface_fourier_radial[p];
            target.inner_surface_fourier_axial[p*2] = cc * sourceCos.inner_surface_fourier_axial[p];
            target.inner_surface_fourier_axial[p*2+1] = ss * sourceSin.inner_surface_fourier_axial[p];
            target.outer_surface_fourier_radial[p*2] = cc * sourceCos.outer_surface_fourier_radial[p];
            target.outer_surface_fourier_radial[p*2+1] = ss * sourceSin.outer_surface_fourier_radial[p];
            target.outer_surface_fourier_axial[p*2] = cc * sourceCos.outer_surface_fourier_axial[p];
            target.outer_surface_fourier_axial[p*2+1] = ss * sourceSin.outer_surface_fourier_axial[p];
        }
    }
    return out;
}

void compareDynamics(const ClapperDynamicsResult& a, const ClapperDynamicsResult& b) {
    require(a.contact_occurred && b.contact_occurred, "rotated clapper case did not contact");
    require(near(a.first_contact_time_seconds, b.first_contact_time_seconds), "first-contact time changed under common rotation");
    require(near(a.separation_time_seconds, b.separation_time_seconds), "separation time changed under common rotation");
    require(near(a.peak_contact_force_n, b.peak_contact_force_n), "peak contact force changed under common rotation");
    require(near(a.first_contact_surface_r_m, b.first_contact_surface_r_m), "contact radius changed under common rotation");
    require(near(a.first_contact_surface_z_m, b.first_contact_surface_z_m), "contact height changed under common rotation");
    require(a.modal_displacement_at_separation.size() == b.modal_displacement_at_separation.size(), "modal displacement size mismatch");
    require(a.modal_velocity_at_separation.size() == b.modal_velocity_at_separation.size(), "modal velocity size mismatch");
    for (std::size_t i=0; i<a.modal_displacement_at_separation.size(); ++i) {
        require(near(a.modal_displacement_at_separation[i], b.modal_displacement_at_separation[i]), "modal release displacement changed under common rotation");
        require(near(a.modal_velocity_at_separation[i], b.modal_velocity_at_separation[i]), "modal release velocity changed under common rotation");
    }
    require(a.contact_force_history_n.size() == b.contact_force_history_n.size(), "contact history size mismatch");
    for (std::size_t i=0; i<a.contact_force_history_n.size(); ++i)
        require(near(a.contact_force_history_n[i], b.contact_force_history_n[i]), "contact force history changed under common rotation");
}
}

int main() {
    try {
        const std::vector<Point2D> inner{{0.30, 0.0}, {0.30, 1.0}};
        const auto baseModes = makeBaseModes();
        constexpr double rotation = 0.731;
        const auto rotatedModes = rotateDegeneratePair(baseModes, rotation);
        ClapperDynamicsConfig baseConfig;
        baseConfig.mass_kg = 0.10;
        baseConfig.moment_of_inertia_kg_m2 = 0.02025;
        baseConfig.centre_of_mass_length_m = 0.45;
        baseConfig.tip_length_m = 0.50;
        baseConfig.tip_radius_m = 0.05;
        baseConfig.initial_angular_velocity_rad_s = 3.6;
        baseConfig.hertz_stiffness_n_per_m_1_5 = 1.0e7;
        baseConfig.integration_step_seconds = 2.0e-6;
        baseConfig.maximum_duration_seconds = 0.5;
        auto rotatedConfig = baseConfig;
        rotatedConfig.strike_azimuth_rad = rotation;

        ClapperDynamicsResult baseResult, rotatedResult;
        require(simulateInternalClapperContact(baseConfig, inner, baseModes, baseResult), "base clapper simulation failed");
        require(simulateInternalClapperContact(rotatedConfig, inner, rotatedModes, rotatedResult), "rotated clapper simulation failed");
        compareDynamics(baseResult, rotatedResult);

        // The acoustic modal fields must obey the same covariance: rotating
        // the modal basis and listener together leaves complex pressure intact.
        constexpr double polar = 0.82, azimuth = -0.37;
        const auto baseField0 = buildDirectionalRadiationField(baseModes[0], inner, polar, 2.0);
        const auto baseField1 = buildDirectionalRadiationField(baseModes[1], inner, polar, 2.0);
        const auto rotatedField0 = buildDirectionalRadiationField(rotatedModes[0], inner, polar, 2.0);
        const auto rotatedField1 = buildDirectionalRadiationField(rotatedModes[1], inner, polar, 2.0);
        const auto q0 = baseResult.modal_velocity_at_separation[0];
        const auto q1 = baseResult.modal_velocity_at_separation[1];
        const auto qr0 = rotatedResult.modal_velocity_at_separation[0];
        const auto qr1 = rotatedResult.modal_velocity_at_separation[1];
        const auto basePressure = q0*baseField0.evaluate(azimuth) + q1*baseField1.evaluate(azimuth);
        const auto rotatedPressure = qr0*rotatedField0.evaluate(azimuth+rotation) +
                                     qr1*rotatedField1.evaluate(azimuth+rotation);
        require(nearComplex(basePressure, rotatedPressure), "complex stereo radiation changed under common clapper/modal/listener rotation");

        std::cout << "clapper_stereo_frame_covariance_contract: PASS (contact trajectory, release modal state, complex radiation)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "clapper_stereo_frame_covariance_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
