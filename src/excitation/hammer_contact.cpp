#include "hammer_contact.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;
using State = std::vector<double>; // striker x,v followed by q,v pairs

struct System {
    const HammerContactConfig& config;
    const std::vector<ContactModalMode>& modes;
    std::vector<double> omega;

    double contactDisplacement(const State& s) const {
        double y = 0.0;
        for (std::size_t i = 0; i < modes.size(); ++i)
            y += modes[i].mode_shape_at_contact * s[2 + 2 * i];
        return y;
    }
    double force(const State& s) const {
        const double delta = s[0] - contactDisplacement(s);
        return delta > 0.0 ? config.hertz_stiffness_n_per_m_1_5 * delta * std::sqrt(delta) : 0.0;
    }
    State derivative(const State& s) const {
        State d(s.size(), 0.0);
        const double f = force(s);
        d[0] = s[1];
        d[1] = -f / config.striker_mass_kg;
        for (std::size_t i = 0; i < modes.size(); ++i) {
            const std::size_t qi = 2 + 2 * i;
            d[qi] = s[qi + 1];
            d[qi + 1] = (modes[i].mode_shape_at_contact * f -
                         modes[i].modal_mass * omega[i] * omega[i] * s[qi]) /
                        modes[i].modal_mass;
        }
        return d;
    }
    double energy(const State& s) const {
        double e = 0.5 * config.striker_mass_kg * s[1] * s[1];
        for (std::size_t i = 0; i < modes.size(); ++i) {
            const std::size_t qi = 2 + 2 * i;
            const double w = omega[i];
            e += 0.5 * modes[i].modal_mass *
                 (s[qi + 1] * s[qi + 1] + w * w * s[qi] * s[qi]);
        }
        const double delta = s[0] - contactDisplacement(s);
        if (delta > 0.0)
            e += 0.4 * config.hertz_stiffness_n_per_m_1_5 * delta * delta * std::sqrt(delta);
        return e;
    }
};

State addScaled(const State& a, const State& b, double scale) {
    State out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) out[i] = a[i] + scale * b[i];
    return out;
}
void rk4Step(const System& sys, State& s, double h) {
    const State k1 = sys.derivative(s);
    const State k2 = sys.derivative(addScaled(s, k1, 0.5 * h));
    const State k3 = sys.derivative(addScaled(s, k2, 0.5 * h));
    const State k4 = sys.derivative(addScaled(s, k3, h));
    for (std::size_t i = 0; i < s.size(); ++i)
        s[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}
} // namespace

bool simulateHertzHammerContact(const HammerContactConfig& c,
                                const std::vector<ContactModalMode>& modes,
                                HammerContactResult& result) {
    result = {};
    if (!std::isfinite(c.striker_mass_kg) || !(c.striker_mass_kg > 0.0) ||
        !std::isfinite(c.initial_velocity_m_per_s) || !(c.initial_velocity_m_per_s > 0.0) ||
        !std::isfinite(c.hertz_stiffness_n_per_m_1_5) || !(c.hertz_stiffness_n_per_m_1_5 > 0.0) ||
        !std::isfinite(c.integration_step_seconds) || !(c.integration_step_seconds > 0.0) ||
        !std::isfinite(c.duration_seconds) || !(c.duration_seconds > 0.0) || modes.empty())
        return false;

    const double stepsReal = std::ceil(c.duration_seconds / c.integration_step_seconds);
    if (!std::isfinite(stepsReal) || stepsReal < 1.0 || stepsReal > 2000000.0) return false;
    System sys{c, modes, {}};
    sys.omega.reserve(modes.size());
    for (const auto& mode : modes) {
        if (!std::isfinite(mode.frequency_hz) || !(mode.frequency_hz > 0.0) ||
            !std::isfinite(mode.modal_mass) || !(mode.modal_mass > 0.0) ||
            !std::isfinite(mode.mode_shape_at_contact)) return false;
        sys.omega.push_back(2.0 * PI * mode.frequency_hz);
    }

    const std::size_t steps = static_cast<std::size_t>(stepsReal);
    const double h = c.duration_seconds / static_cast<double>(steps);
    State s(2 + 2 * modes.size(), 0.0);
    s[1] = c.initial_velocity_m_per_s;
    const double e0 = sys.energy(s);
    if (!(e0 > 0.0) || !std::isfinite(e0)) return false;

    result.initial_energy_joules = e0;
    result.striker_displacement_m.reserve(steps + 1);
    result.striker_velocity_m_per_s.reserve(steps + 1);
    result.contact_force_n.reserve(steps + 1);
    result.total_energy_joules.reserve(steps + 1);
    result.modal_displacement.assign(modes.size(), {});
    result.modal_velocity.assign(modes.size(), {});
    for (auto& v : result.modal_displacement) v.reserve(steps + 1);
    for (auto& v : result.modal_velocity) v.reserve(steps + 1);

    for (std::size_t n = 0; n <= steps; ++n) {
        const double f = sys.force(s);
        const double e = sys.energy(s);
        if (!std::isfinite(f) || !std::isfinite(e) || !std::isfinite(s[0]) || !std::isfinite(s[1])) {
            result = {};
            return false;
        }
        result.striker_displacement_m.push_back(s[0]);
        result.striker_velocity_m_per_s.push_back(s[1]);
        result.contact_force_n.push_back(f);
        result.total_energy_joules.push_back(e);
        if (f > 0.0) result.contact_duration_seconds += h;
        for (std::size_t i = 0; i < modes.size(); ++i) {
            const double q = s[2 + 2 * i];
            if (!std::isfinite(q)) { result = {}; return false; }
            const double qv = s[2 + 2 * i + 1];
            if (!std::isfinite(q) || !std::isfinite(qv)) { result = {}; return false; }
            result.modal_displacement[i].push_back(q);
            result.modal_velocity[i].push_back(qv);
        }
        result.maximum_relative_energy_error = std::max(
            result.maximum_relative_energy_error, std::abs(e - e0) / e0);
        if (n < steps) rk4Step(sys, s, h);
    }
    result.final_energy_joules = result.total_energy_joules.back();
    return true;
}

} // namespace BellOptimizer::Excitation
