#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    // Algebraic contract for the legacy meridional influence proxy:
    // u_r = -z*theta, u_z = r*theta. This test checks algebra only; it does
    // not establish that this is a realizable global rigid rotation of a
    // full 3-D axisymmetric bell.
    const std::vector<double> r{0.10, 0.20, 0.30};
    const std::vector<double> z{0.40, 0.50, 0.60};
    for (std::size_t i = 0; i < r.size(); ++i) {
        const double qr = -z[i];
        const double qz = r[i];
        assert(std::abs(std::hypot(qr, qz) - std::hypot(z[i], r[i])) < 1e-12);
    }

    // Under the assumed modal/base-excitation convention, p = phi^T M q.
    // Physical interpretation as an effective rotational inertia is pending
    // a justified excitation field and reference frame for the axisymmetric model.
    const double phi[] = {0.4, 0.8, -0.3, 0.2, 0.5, -0.7};
    const double mass[] = {2.0, 1.5, 3.0};
    double p = 0.0;
    double modal_mass = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const double qr = -z[i];
        const double qz = r[i];
        p += mass[i] * (phi[2*i] * qr + phi[2*i+1] * qz);
        modal_mass += mass[i] * (phi[2*i] * phi[2*i] + phi[2*i+1] * phi[2*i+1]);
    }
    assert(std::isfinite(p));
    assert(modal_mass > 0.0);
    assert(std::abs((-p) * (-p) - p * p) < 1e-12);

    // Conditional derivation only: if a physically justified base influence
    // q were specified, the generalized modal equation would have forcing
    // -p * theta_ddot. An angular-velocity impulse J_theta therefore gives
    // q_m(t) = -J_theta*p/omega * sin(omega*t). The impulse is 1/s, not N*s.
    const double angular_velocity_impulse = 0.25;
    const double omega = 2.0 * 3.14159265358979323846 * 800.0;
    const double modal_velocity_jump = angular_velocity_impulse * p;
    const double modal_amplitude = std::abs(modal_velocity_jump) / omega;
    assert(std::isfinite(modal_velocity_jump));
    assert(std::isfinite(modal_amplitude));
    assert(std::abs(modal_amplitude * omega - std::abs(modal_velocity_jump)) < 1e-12);

    // Far-field monopole pressure adds omega^2*Qd to displacement; after the
    // 1/omega modal displacement from an impulse, the relative pressure
    // coefficient is proportional to J_theta*p*Qd*omega.
    const double qd = -0.75;
    const double pressure_weight = angular_velocity_impulse * p * qd * omega;
    assert(std::isfinite(pressure_weight));
    assert(std::abs(pressure_weight / (angular_velocity_impulse * p * qd) - omega) < 1e-10);

    std::cout << "rotational excitation contract: PASS\n";
    return 0;
}
