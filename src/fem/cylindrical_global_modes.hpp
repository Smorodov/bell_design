#pragma once

#include "cylindrical_global_assembly.hpp"
#include "sparse_eigensolver.hpp"
#include <vector>

namespace BellOptimizer::FEM {

struct CylindricalConstraintTransform {
    // For each full-system DOF: reduced index, or -1 when constrained to zero.
    // The full DOF value is scale[full] * reduced[reduced_index[full]].
    // Multiple full DOFs may map to one reduced DOF (coupled constraints).
    std::vector<int> reduced_index;
    std::vector<double> scale;
    int reduced_dof_count = 0;
};

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

// Construct exact smooth-vector-field regularity constraints at r=0 for the
// real Fourier basis [0, +1 cosine, -1 sine, +2 cosine, -2 sine, ...].
// The m=1 transverse traces are coupled; no penalty terms are introduced.
CylindricalConstraintTransform make_cylindrical_axis_regularity_transform(
    const BellMeshData& mesh, std::size_t dofs_per_node,
    const std::vector<int>& fourier_basis_labels, double axis_tolerance = 1e-9);

// Solve after exact linear constraint elimination u=Pq. Matrices are reduced
// as P^T K P and P^T M P; returned mode shapes are expanded to full DOF order.
CylindricalModalResult solve_cylindrical_global_modes(
    const CylindricalGlobalSystem& system, const CylindricalConstraintTransform& transform,
    int requested_modes, int max_iterations = 160, double tolerance = 1e-8,
    bool use_rcm_ordering = true,
    const std::vector<std::vector<double>>& initial_modes = {});

} // namespace BellOptimizer::FEM
