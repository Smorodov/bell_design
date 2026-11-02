#ifdef NDEBUG
#undef NDEBUG
#endif

#include "excitation/hammer_contact.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace BellOptimizer::Excitation;

int main() {
    const std::vector<ContactModalMode> modes{
        {420.0, 0.45, 0.8},
        {730.0, 0.30, 0.25}
    };
    HammerContactConfig c;
    c.striker_mass_kg = 0.045;
    c.initial_velocity_m_per_s = 1.2;
    c.hertz_stiffness_n_per_m_1_5 = 2.0e7;
    c.integration_step_seconds = 2.0e-6;
    c.duration_seconds = 0.025;

    HammerContactResult coarse;
    const bool coarseOk = simulateHertzHammerContact(c, modes, coarse);
    assert(coarseOk);
    if (!coarseOk) return 1;
    assert(coarse.striker_displacement_m.size() > 1000);
    assert(coarse.modal_displacement.size() == modes.size());
    assert(coarse.modal_velocity.size() == modes.size());
    assert(coarse.modal_velocity[0].size() == coarse.contact_force_n.size());
    assert(coarse.contact_duration_seconds > 0.0);
    assert(coarse.contact_duration_seconds < c.duration_seconds);
    double peak = 0.0;
    for (double f : coarse.contact_force_n) {
        assert(std::isfinite(f) && f >= 0.0);
        peak = std::max(peak, f);
    }
    assert(peak > 0.0);
    assert(coarse.maximum_relative_energy_error < 2.0e-3);
    assert(std::abs(coarse.final_energy_joules - coarse.initial_energy_joules) /
           coarse.initial_energy_joules < 2.0e-3);

    // Refining the integration step must improve the conservative energy balance.
    c.integration_step_seconds *= 0.5;
    HammerContactResult fine;
    const bool fineOk = simulateHertzHammerContact(c, modes, fine);
    assert(fineOk);
    if (!fineOk) return 1;
    assert(fine.maximum_relative_energy_error < coarse.maximum_relative_energy_error);
    assert(fine.contact_force_n.size() > coarse.contact_force_n.size());
    const double dt = c.duration_seconds / static_cast<double>(fine.contact_force_n.size() - 1);
    double contactImpulse = 0.0;
    for (std::size_t i = 1; i < fine.contact_force_n.size(); ++i)
        contactImpulse += 0.5 * (fine.contact_force_n[i - 1] + fine.contact_force_n[i]) * dt;
    const double strikerMomentumChange = c.striker_mass_kg *
        (fine.striker_velocity_m_per_s.front() - fine.striker_velocity_m_per_s.back());
    assert(std::abs(contactImpulse - strikerMomentumChange) /
           std::max(1.0e-12, std::abs(contactImpulse)) < 2.0e-4);

    HammerContactResult invalid;
    c.striker_mass_kg = 0.0;
    const bool rejectsBadMass = simulateHertzHammerContact(c, modes, invalid);
    assert(!rejectsBadMass);
    if (rejectsBadMass) return 1;
    c.striker_mass_kg = 0.045;
    const bool rejectsNoModes = simulateHertzHammerContact(c, {}, invalid);
    assert(!rejectsNoModes);
    if (rejectsNoModes) return 1;

    std::cout << "hammer_contact_contract: PASS; contact="
              << fine.contact_duration_seconds << " s; peak_force=";
    double finePeak = 0.0;
    for (double f : fine.contact_force_n) finePeak = std::max(finePeak, f);
    std::cout << finePeak << " N; energy_error="
              << fine.maximum_relative_energy_error << "\n";
}
