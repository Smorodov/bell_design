#include "clapper_dynamics.hpp"
#include "internal_clapper_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {
namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kGravity = 9.80665;

struct ModeContact {
    double omega = 0.0;
    double modal_mass = 0.0;
    const ModalMode* mode = nullptr;
};
using State = std::vector<double>; // theta, angular velocity, then q/v pairs

struct ContactEvaluation {
    double force = 0.0;
    double penetration = 0.0;
    Point2D surface_point{};
    double generalized_torque = 0.0;
    std::vector<double> normal_shapes;
};

bool finiteMode(const ModalMode& m, std::size_t profileCount) {
    const std::size_t basis = m.inner_surface_harmonics.size();
    return m.frequency_hz > 0.0 && std::isfinite(m.frequency_hz) &&
           m.modal_mass > 0.0 && std::isfinite(m.modal_mass) && basis > 0 &&
           m.inner_surface_fourier_radial.size() == profileCount * basis &&
           m.inner_surface_fourier_axial.size() == profileCount * basis;
}

double fourierAt(const std::vector<double>& values, std::size_t point,
                 std::size_t pointCount, const std::vector<int>& labels,
                 double azimuth) {
    if (values.size() != pointCount * labels.size()) return 0.0;
    double result = 0.0;
    for (std::size_t b = 0; b < labels.size(); ++b) {
        const int h = labels[b];
        const double v = values[point * labels.size() + b];
        if (h == 0) result += v;
        else if (h > 0) result += v * std::cos(static_cast<double>(h) * azimuth);
        else result += v * std::sin(static_cast<double>(-h) * azimuth);
    }
    return result;
}

struct System {
    const ClapperDynamicsConfig& c;
    const std::vector<Point2D>& profile;
    const std::vector<ModeContact>& modes;

    double tipR(double theta) const { return c.pivot_r_m + c.tip_length_m * std::sin(theta); }
    double tipZ(double theta) const { return c.pivot_z_m + c.tip_length_m * std::cos(theta); }

    double modalNormalDisplacement(const State& s, const SphericalTipContactGeometry& geom,
                                   std::vector<double>* shapes) const {
        double u = 0.0;
        if (shapes) shapes->assign(modes.size(), 0.0);
        const std::size_t i = geom.segment_index;
        const double t = geom.segment_parameter;
        for (std::size_t k = 0; k < modes.size(); ++k) {
            const ModalMode& m = *modes[k].mode;
            const auto& labels = m.inner_surface_harmonics;
            const double ur0 = fourierAt(m.inner_surface_fourier_radial, i, profile.size(), labels, c.strike_azimuth_rad);
            const double ur1 = fourierAt(m.inner_surface_fourier_radial, i + 1, profile.size(), labels, c.strike_azimuth_rad);
            const double uz0 = fourierAt(m.inner_surface_fourier_axial, i, profile.size(), labels, c.strike_azimuth_rad);
            const double uz1 = fourierAt(m.inner_surface_fourier_axial, i + 1, profile.size(), labels, c.strike_azimuth_rad);
            const double ur = (1.0 - t) * ur0 + t * ur1;
            const double uz = (1.0 - t) * uz0 + t * uz1;
            const double shape = ur * geom.outward_normal.r + uz * geom.outward_normal.z;
            if (shapes) (*shapes)[k] = shape;
            u += shape * s[2 + 2 * k];
        }
        return u;
    }

    ContactEvaluation contact(const State& s) const {
        ContactEvaluation out;
        out.normal_shapes.assign(modes.size(), 0.0);
        SphericalTipContactGeometry geom;
        if (!find_spherical_tip_contact_geometry(profile, tipR(s[0]), tipZ(s[0]), c.tip_radius_m, geom))
            return out;
        out.normal_shapes.resize(modes.size());
        const double u = modalNormalDisplacement(s, geom, &out.normal_shapes);
        out.surface_point = geom.surface_point;
        out.penetration = std::max(0.0, u - geom.gap);
        if (!(out.penetration > 0.0)) return out;
        out.force = c.hertz_stiffness_n_per_m_1_5 * out.penetration * std::sqrt(out.penetration);
        // Sphere receives +n*F; generalized torque is F*n dot d(position)/d(theta).
        const double dR = c.tip_length_m * std::cos(s[0]);
        const double dZ = -c.tip_length_m * std::sin(s[0]);
        out.generalized_torque = out.force *
            (geom.outward_normal.r * dR + geom.outward_normal.z * dZ);
        return out;
    }

    State derivative(const State& s) const {
        State d(s.size(), 0.0);
        const ContactEvaluation ce = contact(s);
        d[0] = s[1];
        d[1] = (-c.mass_kg * kGravity * c.centre_of_mass_length_m * std::sin(s[0]) +
                ce.generalized_torque) / c.moment_of_inertia_kg_m2;
        for (std::size_t k = 0; k < modes.size(); ++k) {
            const std::size_t qi = 2 + 2 * k;
            d[qi] = s[qi + 1];
            d[qi + 1] = (-ce.normal_shapes[k] * ce.force -
                         modes[k].modal_mass * modes[k].omega * modes[k].omega * s[qi]) /
                        modes[k].modal_mass;
        }
        return d;
    }

    double energy(const State& s) const {
        double e = 0.5 * c.moment_of_inertia_kg_m2 * s[1] * s[1] +
                   c.mass_kg * kGravity * c.centre_of_mass_length_m * (1.0 - std::cos(s[0]));
        for (std::size_t k = 0; k < modes.size(); ++k) {
            const std::size_t qi = 2 + 2 * k;
            const double w = modes[k].omega;
            e += 0.5 * modes[k].modal_mass *
                 (s[qi + 1] * s[qi + 1] + w * w * s[qi] * s[qi]);
        }
        const ContactEvaluation ce = contact(s);
        if (ce.penetration > 0.0)
            e += 0.4 * c.hertz_stiffness_n_per_m_1_5 *
                 ce.penetration * ce.penetration * std::sqrt(ce.penetration);
        return e;
    }
};

State addScaled(const State& a, const State& b, double scale) {
    State out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) out[i] = a[i] + scale * b[i];
    return out;
}
void rk4(System const& sys, State& s, double h) {
    const State k1 = sys.derivative(s);
    const State k2 = sys.derivative(addScaled(s, k1, 0.5 * h));
    const State k3 = sys.derivative(addScaled(s, k2, 0.5 * h));
    const State k4 = sys.derivative(addScaled(s, k3, h));
    for (std::size_t i = 0; i < s.size(); ++i)
        s[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
}
 } // namespace

bool makeScaledClapperDynamicsConfig(
    double final_bell_height_m, double inner_top_z_m,
    const ClapperScalingParams& p, const BellMaterial& bell_material,
    ClapperDynamicsConfig& config, ClapperScalingResult* summary) {
    config = {};
    if (!(std::isfinite(final_bell_height_m) && final_bell_height_m > 0.0) ||
        !std::isfinite(inner_top_z_m) ||
        !(std::isfinite(p.reference_bell_height_m) && p.reference_bell_height_m > 0.0) ||
        !(std::isfinite(p.reference_mass_kg) && p.reference_mass_kg > 0.0) ||
        !(std::isfinite(p.reference_moment_of_inertia_kg_m2) && p.reference_moment_of_inertia_kg_m2 > 0.0) ||
        !(std::isfinite(p.reference_com_length_m) && p.reference_com_length_m > 0.0) ||
        !(std::isfinite(p.reference_tip_length_m) && p.reference_tip_length_m > 0.0) ||
        !(std::isfinite(p.reference_tip_radius_m) && p.reference_tip_radius_m > 0.0) ||
        !(std::isfinite(p.reference_pivot_below_top_m) && p.reference_pivot_below_top_m >= 0.0) ||
        !std::isfinite(p.reference_initial_angular_velocity_rad_s) ||
        !std::isfinite(p.initial_angle_rad) ||
        !(std::isfinite(p.reference_integration_step_s) && p.reference_integration_step_s > 0.0) ||
        !(std::isfinite(p.reference_maximum_duration_s) && p.reference_maximum_duration_s > 0.0) ||
        !(std::isfinite(p.tip_young_modulus_pa) && p.tip_young_modulus_pa > 0.0) ||
        !std::isfinite(p.tip_poisson) || p.tip_poisson <= -1.0 || p.tip_poisson >= 0.5 ||
        !std::isfinite(p.strike_azimuth_rad) ||
        !(std::isfinite(bell_material.E) && bell_material.E > 0.0) ||
        !std::isfinite(bell_material.nu) || bell_material.nu <= -1.0 || bell_material.nu >= 0.5)
        return false;

    const double s = final_bell_height_m / p.reference_bell_height_m;
    if (!(std::isfinite(s) && s >= 0.05 && s <= 4.0)) return false;
    const double mass = p.reference_mass_kg * s * s * s;
    const double com = p.reference_com_length_m * s;
    const double inertia = p.reference_moment_of_inertia_kg_m2 * s * s * s * s * s;
    const double tip_length = p.reference_tip_length_m * s;
    const double tip_radius = p.reference_tip_radius_m * s;
    const double bell_young_pa = bell_material.E * 1.0e6;
    const double effective_modulus = 1.0 / (
        (1.0 - p.tip_poisson * p.tip_poisson) / p.tip_young_modulus_pa +
        (1.0 - bell_material.nu * bell_material.nu) / bell_young_pa);
    const double hertz = (4.0 / 3.0) * effective_modulus * std::sqrt(tip_radius);
    const double omega = p.reference_initial_angular_velocity_rad_s / std::sqrt(s);
    // Preserve the dimensionless integration resolution under Froude scaling:
    // characteristic impact time scales as sqrt(length), so dt must do the same.
    const double step = p.reference_integration_step_s * std::sqrt(s);
    const double duration = p.reference_maximum_duration_s * std::sqrt(s);
    if (!(std::isfinite(mass) && mass > 0.0 && std::isfinite(inertia) && inertia > 0.0 &&
          std::isfinite(hertz) && hertz > 0.0 && std::isfinite(step) && step > 0.0 &&
          std::isfinite(duration) && duration > 0.0 && duration / step <= 2000000.0)) return false;

    config.mass_kg = mass;
    config.moment_of_inertia_kg_m2 = inertia;
    config.centre_of_mass_length_m = com;
    config.tip_length_m = tip_length;
    config.tip_radius_m = tip_radius;
    config.pivot_r_m = 0.0;
    config.pivot_z_m = inner_top_z_m - p.reference_pivot_below_top_m * s;
    config.initial_angle_rad = p.initial_angle_rad;
    config.initial_angular_velocity_rad_s = omega;
    config.hertz_stiffness_n_per_m_1_5 = hertz;
    config.strike_azimuth_rad = p.strike_azimuth_rad;
    config.integration_step_seconds = step;
    config.maximum_duration_seconds = duration;
    if (summary) {
        summary->linear_scale = s;
        summary->mass_kg = mass;
        summary->moment_of_inertia_kg_m2 = inertia;
        summary->centre_of_mass_length_m = com;
        summary->tip_length_m = tip_length;
        summary->tip_radius_m = tip_radius;
        summary->hertz_stiffness_n_per_m_1_5 = hertz;
    }
    return true;
}

bool simulateInternalClapperContact(const ClapperDynamicsConfig& c,
                                    const std::vector<Point2D>& profile,
                                    const std::vector<ModalMode>& sourceModes,
                                    ClapperDynamicsResult& result) {
    result = {};
    if (profile.size() < 2 || sourceModes.empty() ||
        !std::isfinite(c.mass_kg) || c.mass_kg <= 0.0 ||
        !std::isfinite(c.moment_of_inertia_kg_m2) || c.moment_of_inertia_kg_m2 <= 0.0 ||
        !std::isfinite(c.centre_of_mass_length_m) || c.centre_of_mass_length_m <= 0.0 ||
        !std::isfinite(c.tip_length_m) || c.tip_length_m <= 0.0 ||
        !std::isfinite(c.tip_radius_m) || c.tip_radius_m <= 0.0 ||
        !std::isfinite(c.pivot_r_m) || c.pivot_r_m < 0.0 || !std::isfinite(c.pivot_z_m) ||
        !std::isfinite(c.initial_angle_rad) || std::abs(c.initial_angle_rad) >= 1.2 ||
        !std::isfinite(c.initial_angular_velocity_rad_s) ||
        !std::isfinite(c.hertz_stiffness_n_per_m_1_5) || c.hertz_stiffness_n_per_m_1_5 <= 0.0 ||
        !std::isfinite(c.integration_step_seconds) || c.integration_step_seconds <= 0.0 ||
        !std::isfinite(c.maximum_duration_seconds) || c.maximum_duration_seconds <= 0.0 ||
        !std::isfinite(c.strike_azimuth_rad)) return false;

    const double stepsReal = std::ceil(c.maximum_duration_seconds / c.integration_step_seconds);
    if (!std::isfinite(stepsReal) || stepsReal < 2.0 || stepsReal > 2000000.0) return false;
    std::vector<ModeContact> modes;
    modes.reserve(sourceModes.size());
    for (const auto& m : sourceModes) {
        if (!finiteMode(m, profile.size())) continue;
        modes.push_back({2.0 * kPi * m.frequency_hz, m.modal_mass, &m});
    }
    if (modes.size() < 3) return false;
    System sys{c, profile, modes};
    State s(2 + 2 * modes.size(), 0.0);
    s[0] = c.initial_angle_rad;
    s[1] = c.initial_angular_velocity_rad_s;
    const double e0 = sys.energy(s);
    if (!(e0 > 0.0) || !std::isfinite(e0)) return false;
    result.initial_energy_joules = e0;

    const std::size_t steps = static_cast<std::size_t>(stepsReal);
    const double h = c.maximum_duration_seconds / static_cast<double>(steps);
    bool hadContact = false;
    bool separated = false;
    double maxEnergyError = 0.0;
    for (std::size_t n = 0; n <= steps; ++n) {
        const ContactEvaluation ce = sys.contact(s);
        const double time = static_cast<double>(n) * h;
        result.contact_force_history_n.push_back(ce.force);
        result.angle_history_rad.push_back(s[0]);
        result.peak_contact_force_n = std::max(result.peak_contact_force_n, ce.force);
        if (ce.force > 1.0e-8) {
            if (!hadContact) {
                hadContact = true;
                result.first_contact_time_seconds = time;
                result.first_contact_surface_r_m = ce.surface_point.r;
                result.first_contact_surface_z_m = ce.surface_point.z;
            }
        } else if (hadContact && !separated) {
            separated = true;
            result.separation_time_seconds = time;
            result.separation_angle_rad = s[0];
            result.separation_angular_velocity_rad_s = s[1];
            result.modal_displacement_at_separation.resize(modes.size());
            result.modal_velocity_at_separation.resize(modes.size());
            for (std::size_t k = 0; k < modes.size(); ++k) {
                result.modal_displacement_at_separation[k] = s[2 + 2 * k];
                result.modal_velocity_at_separation[k] = s[3 + 2 * k];
            }
            result.final_energy_joules = sys.energy(s);
            result.modal_energy_at_separation_joules.resize(modes.size(), 0.0);
            for (std::size_t k = 0; k < modes.size(); ++k) {
                const std::size_t qi = 2 + 2 * k;
                const double velocity = s[qi + 1];
                const double displacement = s[qi];
                const double omega = modes[k].omega;
                const double modalEnergy = 0.5 * modes[k].modal_mass *
                    (velocity * velocity + omega * omega * displacement * displacement);
                result.modal_energy_at_separation_joules[k] = modalEnergy;
                result.bell_modal_energy_joules += modalEnergy;
            }
            result.clapper_energy_at_separation_joules =
                0.5 * c.moment_of_inertia_kg_m2 * s[1] * s[1] +
                c.mass_kg * kGravity * c.centre_of_mass_length_m * (1.0 - std::cos(s[0]));
            result.bell_energy_transfer_fraction = result.initial_energy_joules > 0.0
                ? result.bell_modal_energy_joules / result.initial_energy_joules : 0.0;
            result.contact_duration_seconds = result.separation_time_seconds - result.first_contact_time_seconds;
            result.maximum_relative_energy_error = maxEnergyError;
            result.contact_occurred = true;
            return true;
        }
        const double energy = sys.energy(s);
        if (std::isfinite(energy)) maxEnergyError = std::max(maxEnergyError, std::abs(energy - e0) / e0);
        if (n < steps) rk4(sys, s, h);
    }
    return false; // no complete first-contact/separation event within the configured horizon
}

} // namespace BellOptimizer::Excitation
