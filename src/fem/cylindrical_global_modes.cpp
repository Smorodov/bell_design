#include "cylindrical_global_modes.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::FEM {
namespace {
detail::SparseMatrix reduce_matrix(const detail::SparseMatrix& source,
                                   const std::vector<int>& free_dofs,
                                   const std::vector<int>& old_to_new) {
    std::vector<detail::SparseTriplet> entries;
    entries.reserve(source.values.size());
    for (std::size_t i = 0; i < free_dofs.size(); ++i) {
        const int old_row = free_dofs[i];
        for (std::size_t p = source.row_offsets[static_cast<std::size_t>(old_row)];
             p < source.row_offsets[static_cast<std::size_t>(old_row) + 1u]; ++p) {
            const int new_col = old_to_new[static_cast<std::size_t>(source.columns[p])];
            if (new_col >= 0) entries.push_back({static_cast<int>(i), new_col, source.values[p]});
        }
    }
    return detail::SparseMatrix::from_triplets(static_cast<int>(free_dofs.size()), std::move(entries));
}
}

CylindricalModalResult solve_cylindrical_global_modes(
    const CylindricalGlobalSystem& system, const std::vector<int>& fixed_dofs,
    int requested_modes, int max_iterations, double tolerance, bool use_rcm_ordering,
    const std::vector<std::vector<double>>& initial_modes) {
    if (system.stiffness.n <= 0 || system.stiffness.n != system.mass.n)
        throw std::invalid_argument("cylindrical modes: inconsistent global matrices");
    if (requested_modes <= 0 || max_iterations <= 0 || !(tolerance > 0.0) || !std::isfinite(tolerance))
        throw std::invalid_argument("cylindrical modes: invalid solver controls");
    const int n = system.stiffness.n;
    std::vector<bool> fixed(static_cast<std::size_t>(n), false);
    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= n) throw std::out_of_range("cylindrical modes: fixed DOF out of range");
        if (fixed[static_cast<std::size_t>(dof)]) throw std::invalid_argument("cylindrical modes: duplicate fixed DOF");
        fixed[static_cast<std::size_t>(dof)] = true;
    }
    std::vector<int> free_dofs, old_to_new(static_cast<std::size_t>(n), -1);
    free_dofs.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) if (!fixed[static_cast<std::size_t>(i)]) {
        old_to_new[static_cast<std::size_t>(i)] = static_cast<int>(free_dofs.size());
        free_dofs.push_back(i);
    }
    if (free_dofs.size() <= static_cast<std::size_t>(requested_modes))
        throw std::invalid_argument("cylindrical modes: too few unconstrained DOFs");
    const auto k = reduce_matrix(system.stiffness, free_dofs, old_to_new);
    const auto m = reduce_matrix(system.mass, free_dofs, old_to_new);
    std::vector<std::vector<double>> reduced_initial;
    for (const auto& full : initial_modes) {
        if (full.size() != static_cast<std::size_t>(n)) continue;
        std::vector<double> reduced(free_dofs.size(), 0.0);
        for (std::size_t i = 0; i < free_dofs.size(); ++i) reduced[i] = full[static_cast<std::size_t>(free_dofs[i])];
        reduced_initial.push_back(std::move(reduced));
    }
    const auto eig = detail::solve_lowest_sparse_modes(k, m, requested_modes, max_iterations,
                                                        tolerance, use_rcm_ordering, reduced_initial);
    CylindricalModalResult out;
    out.converged = eig.converged;
    out.relative_residual = eig.final_relative_residual;
    out.iterations = eig.iterations;
    if (!eig.converged) return out;
    constexpr double two_pi = 6.283185307179586476925286766559;
    for (std::size_t mode = 0; mode < eig.eigenvalues.size(); ++mode) {
        const double lambda = eig.eigenvalues[mode];
        const double spectral_scale = eig.eigenvalues.empty() ? 0.0 : std::abs(eig.eigenvalues.back());
        const double zero_tolerance = std::max(spectral_scale * 1e-12, std::numeric_limits<double>::min());
        if (!(lambda > zero_tolerance) || !std::isfinite(lambda)) continue;
        std::vector<double> full(static_cast<std::size_t>(n), 0.0);
        for (std::size_t i = 0; i < free_dofs.size(); ++i)
            full[static_cast<std::size_t>(free_dofs[i])] = eig.vectors[mode][i];
        out.eigenvalues.push_back(lambda);
        out.frequencies_hz.push_back(std::sqrt(lambda) / two_pi);
        out.mode_shapes.push_back(std::move(full));
        if (out.frequencies_hz.size() == static_cast<std::size_t>(requested_modes)) break;
    }
    if (out.frequencies_hz.size() < static_cast<std::size_t>(requested_modes)) out.converged = false;
    return out;
}
} // namespace BellOptimizer::FEM
