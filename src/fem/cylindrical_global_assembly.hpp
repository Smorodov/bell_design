#pragma once

#include "core/bell_types.h"
#include "sparse_matrix.hpp"
#include <vector>

namespace BellOptimizer::FEM {

struct CylindricalGlobalSystem {
    detail::SparseMatrix stiffness;
    detail::SparseMatrix mass;
    std::size_t basis_count = 0;
    std::size_t dofs_per_node = 0;
};

// Assemble the cylindrical Q4 Fourier-basis elements into sparse global K/M.
// ellipticity maps the circular reference section to an ellipse (q=b/a) in plan view.
// Global DOF order: (node, displacement component [ur,uphi,uz], Fourier basis).
// Mesh element node indices are one-based mesh node indices (same convention as Element in BellMeshData).
CylindricalGlobalSystem assemble_cylindrical_global_system(
    const BellMeshData& mesh, double young_modulus_pa, double poisson_ratio,
    double density_kg_m3, const std::vector<int>& harmonics = {0, 1, 2},
    unsigned angular_samples = 64, double zero_tolerance = 0.0, double ellipticity = 1.0);

} // namespace BellOptimizer::FEM
