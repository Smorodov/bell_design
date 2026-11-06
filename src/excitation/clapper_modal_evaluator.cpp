#include "clapper_modal_evaluator.hpp"

#include <algorithm>
#include <cmath>

namespace BellOptimizer::Excitation {
namespace {
bool finite(double x) { return std::isfinite(x); }
}

bool evaluateClapperModalState(
    const ClapperModalEvaluationConfig& config,
    const std::vector<Point2D>& profile,
    const std::vector<ModalMode>& modes,
    ClapperModalEvaluationResult& result) {
    result = {};
    const double span = config.mouth_z_m - config.top_z_m;
    if (profile.size() < 2 || modes.size() < 3 || !finite(span) || std::abs(span) < 1.0e-9 ||
        !finite(config.top_z_m) || !finite(config.mouth_z_m) ||
        !finite(config.admissible_strike_fraction_min) ||
        !finite(config.admissible_strike_fraction_max) ||
        config.admissible_strike_fraction_min < 0.0 ||
        config.admissible_strike_fraction_max > 1.0 ||
        config.admissible_strike_fraction_min > config.admissible_strike_fraction_max) return false;

    for (const auto& p : profile) if (!finite(p.r) || !finite(p.z) || p.r < 0.0) return false;
    for (const auto& m : modes) {
        if (!finite(m.frequency_hz) || m.frequency_hz <= 0.0 ||
            !finite(m.modal_mass) || m.modal_mass <= 0.0 ||
            !finite(m.acoustic_pressure_coefficient) ||
            m.inner_surface_harmonics.empty() ||
            m.inner_surface_fourier_radial.size() != profile.size() * m.inner_surface_harmonics.size() ||
            m.inner_surface_fourier_axial.size() != profile.size() * m.inner_surface_harmonics.size()) return false;
    }

    if (!simulateInternalClapperContact(config.clapper, profile, modes, result.contact) ||
        !result.contact.contact_occurred || result.contact.modal_displacement_at_separation.size() != modes.size() ||
        result.contact.modal_velocity_at_separation.size() != modes.size() ||
        result.contact.angle_history_rad.empty()) return false;

    // Use the actual closest point on the inner wall, not the spherical head
    // centre, so the reported strike fraction corresponds to the physical
    // contact location on the bell profile.
    result.first_contact_z_m = result.contact.first_contact_surface_z_m;
    result.strike_height_fraction = (result.first_contact_z_m - config.top_z_m) / span;
    result.strike_height_admissible =
        result.strike_height_fraction >= config.admissible_strike_fraction_min &&
        result.strike_height_fraction <= config.admissible_strike_fraction_max;

    result.modal_pressure_cosine.resize(modes.size());
    result.modal_pressure_sine.resize(modes.size());
    result.modal_pressure_magnitude.resize(modes.size());
    constexpr double two_pi = 6.283185307179586476925286766559;
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const double omega = two_pi * modes[i].frequency_hz;
        const double q = result.contact.modal_displacement_at_separation[i];
        const double v = result.contact.modal_velocity_at_separation[i];
        // Free response q(t)=q0 cos(wt)+(v0/w) sin(wt). The signed
        // radiation coefficient is retained in both quadratures (phase).
        const double c = modes[i].acoustic_pressure_coefficient;
        result.modal_pressure_cosine[i] = c * omega * q;
        result.modal_pressure_sine[i] = c * v;
        result.modal_pressure_magnitude[i] = std::hypot(
            result.modal_pressure_cosine[i], result.modal_pressure_sine[i]);
        if (!finite(result.modal_pressure_magnitude[i])) return false;
    }
    return true;
}

} // namespace BellOptimizer::Excitation
