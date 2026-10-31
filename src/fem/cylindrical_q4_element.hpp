#pragma once

#include "math/cylindrical_strain.hpp"
#include <array>
#include <vector>

namespace BellOptimizer::FEM {

struct MeridionalNode { double r = 0.0; double z = 0.0; };
struct CylindricalQ4ElementMatrices {
    // DOF ordering: ((node * 3 + cylindrical_component) * basis_count + basis_index).
    // Components are [u_r, u_phi, u_z]; basis order is the returned angular basis.
    std::vector<int> harmonics;
    std::vector<bool> sine_basis;
    std::vector<std::vector<double>> stiffness;
    std::vector<std::vector<double>> mass;
};

// Integrates a Q4 meridional element over the revolved 3-D volume using a real
// Fourier basis in phi. This is a local element kernel, not yet a global mesh
// assembler or a production modal solver.
CylindricalQ4ElementMatrices assemble_cylindrical_q4(
    const std::array<MeridionalNode, 4>& nodes,
    double young_modulus_pa, double poisson_ratio, double density_kg_m3,
    const std::vector<int>& even_harmonics = {0, 2},
    unsigned angular_samples = 64, double ellipticity = 1.0);

} // namespace BellOptimizer::FEM
