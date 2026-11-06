#include "excitation/clapper_dynamics.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    try {
        // Vertical inner wall at r=0.30 m; cavity lies toward decreasing r.
        const std::vector<Point2D> inner{{0.30, 0.0}, {0.30, 1.0}};
        std::vector<ModalMode> modes(3);
        const double frequencies[] = {180.0, 270.0, 410.0};
        for (std::size_t i = 0; i < modes.size(); ++i) {
            modes[i].frequency_hz = frequencies[i];
            modes[i].modal_mass = 1.0;
            modes[i].inner_surface_harmonics = {0};
            modes[i].inner_surface_fourier_radial = {0.001, 0.001};
            modes[i].inner_surface_fourier_axial = {0.0, 0.0};
        }
        ClapperDynamicsConfig c;
        c.mass_kg = 0.10;
        c.moment_of_inertia_kg_m2 = 0.02025;
        c.centre_of_mass_length_m = 0.45;
        c.tip_length_m = 0.50;
        c.tip_radius_m = 0.05;
        c.pivot_z_m = 0.0;
        c.initial_angle_rad = 0.0;
        c.initial_angular_velocity_rad_s = 3.6;
        c.hertz_stiffness_n_per_m_1_5 = 1.0e7;
        c.integration_step_seconds = 2.0e-6;
        c.maximum_duration_seconds = 0.5;

        ClapperDynamicsResult result;
        require(simulateInternalClapperContact(c, inner, modes, result),
                "pendulum did not complete a contact and release");
        require(result.contact_occurred, "contact event was not recorded");
        require(result.peak_contact_force_n > 0.0, "contact force is not positive");
        require(std::abs(result.first_contact_surface_r_m - 0.30) < 1.0e-9 &&
                result.first_contact_surface_z_m >= 0.0 && result.first_contact_surface_z_m <= 1.0,
                "first-contact surface coordinates were not recorded on the inner wall");
        require(result.contact_duration_seconds > 0.0 && result.contact_duration_seconds < 0.02,
                "contact duration is outside the expected short-impact range");
        require(result.maximum_relative_energy_error < 2.0e-4,
                "coupled pendulum/contact/modal energy drift is too large");
        require(result.modal_displacement_at_separation.size() == modes.size() &&
                result.modal_velocity_at_separation.size() == modes.size(),
                "modal release state has the wrong size");
        std::cout << "clapper_dynamics_contract: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "clapper_dynamics_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
