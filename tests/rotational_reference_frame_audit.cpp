#include <cmath>
#include <iostream>
#include <vector>

int main() {
    // Discrete analogue of p = phi^T M q with q_r = -(z-z_ref), q_z = r.
    // Translating z_ref by dz changes q_r by dz, so p changes by dz times
    // the radial translational participation. This is expected, not invariant.
    const double r[] = {0.12, 0.20, 0.31};
    const double z[] = {0.25, 0.47, 0.73};
    const double phi_r[] = {0.4, -0.2, 0.7};
    const double phi_z[] = {-0.1, 0.8, 0.3};
    const double mass[] = {1.2, 2.1, 0.9};
    auto participation = [&](double z_ref) {
        double p = 0.0;
        for (int i = 0; i < 3; ++i)
            p += mass[i] * (phi_r[i] * (-(z[i] - z_ref)) + phi_z[i] * r[i]);
        return p;
    };
    double radial_participation = 0.0;
    for (int i = 0; i < 3; ++i) radial_participation += mass[i] * phi_r[i];
    const double dz = 0.37;
    const double p0 = participation(0.0);
    const double p1 = participation(dz);
    const double predicted_delta = dz * radial_participation;
    const double error = std::abs((p1 - p0) - predicted_delta);
    if (error > 1e-12) return 2;
    std::cout << "rotational_reference_frame_audit: PASS\n"
              << "p_at_zref_0=" << p0 << '\n'
              << "p_at_zref_shifted=" << p1 << '\n'
              << "participation_delta=" << (p1 - p0) << '\n'
              << "predicted_delta=" << predicted_delta << '\n'
              << "conclusion=reference_origin_must_be_explicit; participation_is_not_origin_invariant\n";
    return 0;
}
