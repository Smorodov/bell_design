#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// The legacy q=(-z,r) field is audited here as a mathematical influence field.
// The test deliberately does NOT assign it physical 3-D rigid-rotation meaning.
int main() {
    const std::vector<double> mass{2.0, 3.0, 5.0, 7.0};
    const std::vector<double> ur{0.10, -0.20, 0.30, -0.40};
    const std::vector<double> uz{0.50,  0.60, -0.70, 0.80};
    const std::vector<double> r {0.20, 0.40, 0.60, 0.80};
    const std::vector<double> z {1.00, 1.50, 2.00, 2.50};

    auto participation = [&](double z_ref) {
        double p = 0.0;
        for (std::size_t i = 0; i < mass.size(); ++i)
            p += mass[i] * (ur[i] * (-(z[i] - z_ref)) + uz[i] * r[i]);
        return p;
    };

    const double p0 = participation(0.0);
    const double pref = participation(1.25);
    double radial_participation = 0.0;
    for (std::size_t i = 0; i < mass.size(); ++i)
        radial_participation += mass[i] * ur[i];

    // Exact reference-shift identity for q_rz=(-z,r):
    // p(z_ref) = p(0) + z_ref * (phi^T M q_r), q_r=(1,0).
    assert(std::abs(pref - (p0 + 1.25 * radial_participation)) < 1e-12);

    // The dimensions of p^2 are those of a generalized rotational inertia
    // only if q is interpreted as displacement per unit angle. This is a
    // dimensional statement, not proof that the field is a realizable 3-D
    // rigid rotation of the full axisymmetric body.
    const double effective_inertia_proxy = p0 * p0;
    assert(std::isfinite(effective_inertia_proxy));

    // A true axisymmetric rigid translation has an origin-independent field;
    // the q=(-z,r) proxy does not. This is the intended diagnostic distinction.
    assert(std::abs(pref - p0) > 1e-12);

    std::cout << "meridional influence audit: PASS\n";
    std::cout << "p0=" << p0 << "\n";
    std::cout << "p_shifted=" << pref << "\n";
    std::cout << "radial_participation=" << radial_participation << "\n";
    std::cout << "effective_inertia_proxy=" << effective_inertia_proxy << "\n";
    return 0;
}
