#pragma once

#include "cylindrical_global_modes.hpp"
#include <cstddef>
#include <vector>

namespace BellOptimizer::FEM {

// Mass-weighted modal overlap in a common DOF space. The metric matrix is
// supplied explicitly so tracking is independent of eigenvector scaling/sign.
double mass_weighted_modal_overlap(const CylindricalModalResult& a, std::size_t mode_a,
                                   const CylindricalModalResult& b, std::size_t mode_b,
                                   const detail::SparseMatrix& metric_mass);

// Minimum principal overlap between two equal-sized modal subspaces. Values
// near one mean that the subspaces are nearly identical; the result is the
// smallest singular value of the mass-weighted cross-overlap matrix.
double mass_weighted_subspace_similarity(const CylindricalModalResult& a,
                                          const std::vector<std::size_t>& indices_a,
                                          const CylindricalModalResult& b,
                                          const std::vector<std::size_t>& indices_b,
                                          const detail::SparseMatrix& metric_mass);

double mass_weighted_subspace_similarity(const CylindricalModalResult& a,
                                          std::size_t first_a,
                                          const CylindricalModalResult& b,
                                          std::size_t first_b,
                                          std::size_t dimension,
                                          const detail::SparseMatrix& metric_mass);

// Exact embedding of modal vectors from one Fourier basis into a larger nested
// Fourier basis. The meridional mesh and three displacement components must be
// unchanged; missing target harmonics are allowed, but source harmonics must
// all be present in the target. This is the bridge used for harmonic-family
// continuation: it changes coordinates, not the physical modal vector.
std::vector<std::vector<double>> prolongate_cylindrical_modes(
    const CylindricalModalResult& source,
    const std::vector<int>& source_harmonics,
    const std::vector<int>& target_harmonics);

// Principal-angle similarity between two already embedded modal subspaces in
// one common DOF space. Currently supports one- and two-dimensional families.
double mass_weighted_subspace_similarity(
    const std::vector<std::vector<double>>& a,
    const std::vector<std::vector<double>>& b,
    const detail::SparseMatrix& metric_mass);

} // namespace BellOptimizer::FEM
