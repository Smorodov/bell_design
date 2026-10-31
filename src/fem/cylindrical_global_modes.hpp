#pragma once

#include "cylindrical_global_assembly.hpp"
#include "sparse_eigensolver.hpp"
#include <vector>

namespace BellOptimizer::FEM {

struct CylindricalModalResult {
    bool converged = false;
    std::vector<double> eigenvalues;       // lambda = omega^2, in s^-2
    std::vector<double> frequencies_hz;
    std::vector<std::vector<double>> mode_shapes; // full global DOF layout; constrained entries are zero
    double relative_residual = 0.0;
    int iterations = 0;
};

// Solves K x = lambda M x after eliminating fixed global DOFs. Fixed DOFs are
// indices in the original assembled system; no penalty stiffness is added.
CylindricalModalResult solve_cylindrical_global_modes(
    const CylindricalGlobalSystem& system, const std::vector<int>& fixed_dofs,
    int requested_modes, int max_iterations = 160, double tolerance = 1e-8,
    bool use_rcm_ordering = true,
    const std::vector<std::vector<double>>& initial_modes = {});

} // namespace BellOptimizer::FEM
