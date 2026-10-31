#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// This audit deliberately does not add a second excitation mechanism to the
// production model. It establishes the mathematical distinction between the
// current distributed RZ/base-rotation excitation and a localized hammer load.
int main() {
    // A mass-normalized mode is represented here by nodal displacement values.
    // For a localized impulse J_F at a point/node s, the generalized modal
    // impulse is I_m = J_F * phi_s (component-wise dot product with force
    // direction). It depends on the strike location.
    const std::vector<double> phi_r{0.10, -0.20, 0.40};
    const std::vector<double> phi_z{0.30,  0.50, -0.10};
    const std::vector<double> z{0.20, 0.50, 0.80};
    const std::vector<double> r{0.10, 0.20, 0.30};

    const double force_impulse = 1.0;
    const double local_radial = force_impulse * phi_r[1];
    const double local_axial  = force_impulse * phi_z[1];
    assert(std::abs(local_radial + 0.20) < 1e-12);
    assert(std::abs(local_axial - 0.50) < 1e-12);

    // The current RZ model instead uses the distributed kinematic influence
    // q_theta=(-z,r), giving p=phi^T M q_theta. With lumped nodal masses this
    // becomes a weighted sum over the whole structure, not a value at one node.
    const std::vector<double> mass{1.0, 2.0, 3.0};
    double p_rz = 0.0;
    for (std::size_t i = 0; i < mass.size(); ++i)
        p_rz += mass[i] * (phi_r[i] * (-z[i]) + phi_z[i] * r[i]);

    // There is no universal scalar that can turn the localized modal impulse
    // into the RZ participation for all modes: the former samples one spatial
    // point, while the latter is a mass-weighted global projection.
    const double other_strike = force_impulse * phi_r[2];
    assert(std::abs(other_strike - local_radial) > 1e-12);
    assert(std::abs(p_rz) > 1e-12);

    // For a finite hammer contact patch S_h, the correct axisymmetric analogue
    // is a ring/patch load and its modal impulse is an integral of phi·f dS.
    // It approaches the point-load expression as the patch shrinks. This is
    // the proper future extension; it is not equivalent to q_theta.
    const double patch_weight_a = 0.25;
    const double patch_weight_b = 0.75;
    const double patch_modal = patch_weight_a * phi_z[0] + patch_weight_b * phi_z[1];
    assert(std::isfinite(patch_modal));
    assert(std::abs(patch_modal - local_axial) > 1e-12);

    std::cout << "localized-vs-RZ excitation audit: PASS\n";
    std::cout << "rz_participation_proxy=" << p_rz << "\n";
    std::cout << "local_strike_node_1_radial=" << local_radial << "\n";
    std::cout << "local_strike_node_1_axial=" << local_axial << "\n";
    std::cout << "finite_patch_axial_projection=" << patch_modal << "\n";
    return 0;
}
