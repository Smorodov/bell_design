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
    auto eig = detail::solve_lowest_sparse_modes(k, m, requested_modes, max_iterations,
                                                 tolerance, use_rcm_ordering, reduced_initial);
    CylindricalModalResult out;
    out.converged = eig.converged;
    out.relative_residual = eig.final_relative_residual;
    out.iterations = eig.iterations;
    if (!eig.converged) return out;

    // The reduced problem can contain numerically null modes that are valid
    // eigenvectors of the generalized solve but are intentionally removed
    // below by the physical-positive-eigenvalue filter.  Near a circular
    // limit this can consume one requested slot.  Retry with a small spectral
    // surplus rather than returning fewer physical modes than requested.  The
    // normal path is unchanged; the second solve is only paid when filtering
    // actually removes a mode.
    auto count_positive_modes = [&](const detail::SparseEigenResult& candidate) {
        if (candidate.eigenvalues.empty()) return 0;
        const double spectral_scale = std::abs(candidate.eigenvalues.back());
        const double zero_tolerance = std::max(
            spectral_scale * 1e-12, std::numeric_limits<double>::min());
        int count = 0;
        for (double lambda : candidate.eigenvalues) {
            if (lambda > zero_tolerance && std::isfinite(lambda)) ++count;
        }
        return count;
    };
    if (count_positive_modes(eig) < requested_modes) {
        const int retry_modes = requested_modes + std::max(2, requested_modes / 4);
        if (free_dofs.size() > static_cast<std::size_t>(retry_modes)) {
            auto retry = detail::solve_lowest_sparse_modes(
                k, m, retry_modes, max_iterations, tolerance,
                use_rcm_ordering, reduced_initial);
            if (retry.converged && count_positive_modes(retry) >= requested_modes) {
                eig = std::move(retry);
                out.relative_residual = eig.final_relative_residual;
                out.iterations += eig.iterations;
            }
        }
    }
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


CylindricalConstraintTransform make_cylindrical_axis_regularity_transform(
    const BellMeshData& mesh, std::size_t dofs_per_node,
    const std::vector<int>& labels, double axis_tolerance) {
    if (dofs_per_node == 0 || labels.empty() || labels.front() != 0 ||
        !(axis_tolerance >= 0.0) || !std::isfinite(axis_tolerance))
        throw std::invalid_argument("axis regularity: invalid layout or tolerance");
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (labels[i] < 0 || (i && labels[i] <= labels[i - 1u]))
            throw std::invalid_argument("axis regularity: Fourier labels must be increasing and nonnegative");
    }
    std::vector<int> basis_labels;
    basis_labels.reserve(2u * labels.size());
    for (int h : labels) {
        basis_labels.push_back(h);
        if (h > 0) basis_labels.push_back(-h);
    }
    if (basis_labels.size() != labels.size() + static_cast<std::size_t>(
            std::count_if(labels.begin(), labels.end(), [](int h) { return h > 0; })))
        throw std::logic_error("axis regularity: internal Fourier basis mismatch");
    // The assembly's basis count includes both cosine and sine for h>0.
    const std::size_t basis_count = basis_labels.size();
    if (dofs_per_node != 3u * basis_count)
        throw std::invalid_argument("axis regularity: expanded Fourier basis size mismatch");
    if (mesh.nodes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / dofs_per_node)
        throw std::length_error("axis regularity: full DOF count exceeds solver index range");

    const int full_count = static_cast<int>(mesh.nodes.size() * dofs_per_node);
    CylindricalConstraintTransform out;
    out.reduced_index.assign(static_cast<std::size_t>(full_count), -1);
    out.scale.assign(static_cast<std::size_t>(full_count), 0.0);
    auto assign = [&](int full, int reduced, double scale) {
        out.reduced_index[static_cast<std::size_t>(full)] = reduced;
        out.scale[static_cast<std::size_t>(full)] = scale;
    };
    auto find_basis = [&](int label) -> std::size_t {
        const auto it = std::find(basis_labels.begin(), basis_labels.end(), label);
        if (it == basis_labels.end()) throw std::invalid_argument("axis regularity: missing Fourier basis label");
        return static_cast<std::size_t>(it - basis_labels.begin());
    };
    const std::size_t b0 = find_basis(0);
    int next = 0;
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        const bool on_axis = std::abs(mesh.nodes[node].r) < axis_tolerance;
        const int base = static_cast<int>(node * dofs_per_node);
        if (!on_axis) {
            for (std::size_t c = 0; c < dofs_per_node; ++c)
                assign(base + static_cast<int>(c), next++, 1.0);
            continue;
        }
        // m=0: transverse components vanish, axial displacement is free.
        assign(base + static_cast<int>(2u * basis_count + b0), next++, 1.0);
        // m=1: preserve the two constant Cartesian transverse translations.
        // ur(1c)=a, uphi(1s)=-a; ur(1s)=b, uphi(1c)=b.
        const auto cosine_one = std::find(basis_labels.begin(), basis_labels.end(), 1);
        const auto sine_one = std::find(basis_labels.begin(), basis_labels.end(), -1);
        if (cosine_one != basis_labels.end() || sine_one != basis_labels.end()) {
            if (cosine_one == basis_labels.end() || sine_one == basis_labels.end())
                throw std::invalid_argument("axis regularity: incomplete m=1 cosine/sine pair");
            const std::size_t c1 = static_cast<std::size_t>(cosine_one - basis_labels.begin());
            const std::size_t s1 = static_cast<std::size_t>(sine_one - basis_labels.begin());
            const int a = next++;
            const int b = next++;
            assign(base + static_cast<int>(0u * basis_count + c1), a, 1.0);
            assign(base + static_cast<int>(1u * basis_count + s1), a, -1.0);
            assign(base + static_cast<int>(0u * basis_count + s1), b, 1.0);
            assign(base + static_cast<int>(1u * basis_count + c1), b, 1.0);
        }
        // m=1 axial traces and every m>=2 component vanish at the axis.
    }
    out.reduced_dof_count = next;
    return out;
}

CylindricalModalResult solve_cylindrical_global_modes(
    const CylindricalGlobalSystem& system, const CylindricalConstraintTransform& transform,
    int requested_modes, int max_iterations, double tolerance, bool use_rcm_ordering,
    const std::vector<std::vector<double>>& initial_modes) {
    if (system.stiffness.n <= 0 || system.stiffness.n != system.mass.n)
        throw std::invalid_argument("cylindrical modes: inconsistent global matrices");
    if (requested_modes <= 0 || max_iterations <= 0 || !(tolerance > 0.0) || !std::isfinite(tolerance))
        throw std::invalid_argument("cylindrical modes: invalid solver controls");
    const int n = system.stiffness.n;
    if (transform.reduced_index.size() != static_cast<std::size_t>(n) ||
        transform.scale.size() != static_cast<std::size_t>(n) || transform.reduced_dof_count <= requested_modes)
        throw std::invalid_argument("cylindrical modes: invalid constraint transform dimensions");
    std::vector<int> multiplicity(static_cast<std::size_t>(transform.reduced_dof_count), 0);
    for (int i = 0; i < n; ++i) {
        const int reduced = transform.reduced_index[static_cast<std::size_t>(i)];
        const double scale = transform.scale[static_cast<std::size_t>(i)];
        if (reduced < -1 || reduced >= transform.reduced_dof_count || !std::isfinite(scale) ||
            ((reduced < 0) != (scale == 0.0)))
            throw std::invalid_argument("cylindrical modes: malformed constraint transform");
        if (reduced >= 0) ++multiplicity[static_cast<std::size_t>(reduced)];
    }
    if (std::any_of(multiplicity.begin(), multiplicity.end(), [](int count) { return count == 0; }))
        throw std::invalid_argument("cylindrical modes: transform contains an unused reduced DOF");

    auto reduce = [&](const detail::SparseMatrix& source) {
        std::vector<detail::SparseTriplet> entries;
        entries.reserve(source.values.size());
        for (int row = 0; row < n; ++row) {
            const int rr = transform.reduced_index[static_cast<std::size_t>(row)];
            if (rr < 0) continue;
            for (std::size_t p = source.row_offsets[static_cast<std::size_t>(row)];
                 p < source.row_offsets[static_cast<std::size_t>(row) + 1u]; ++p) {
                const int col = source.columns[p];
                const int rc = transform.reduced_index[static_cast<std::size_t>(col)];
                if (rc < 0) continue;
                const double value = source.values[p] * transform.scale[static_cast<std::size_t>(row)] *
                                     transform.scale[static_cast<std::size_t>(col)];
                if (value != 0.0) entries.push_back({rr, rc, value});
            }
        }
        return detail::SparseMatrix::from_triplets(transform.reduced_dof_count, std::move(entries));
    };
    const auto k = reduce(system.stiffness);
    const auto m = reduce(system.mass);

    std::vector<std::vector<double>> reduced_initial;
    for (const auto& full : initial_modes) {
        if (full.size() != static_cast<std::size_t>(n)) continue;
        std::vector<double> reduced(static_cast<std::size_t>(transform.reduced_dof_count), 0.0);
        for (int i = 0; i < n; ++i) {
            const int ri = transform.reduced_index[static_cast<std::size_t>(i)];
            if (ri >= 0) reduced[static_cast<std::size_t>(ri)] +=
                transform.scale[static_cast<std::size_t>(i)] * full[static_cast<std::size_t>(i)] /
                static_cast<double>(multiplicity[static_cast<std::size_t>(ri)]);
        }
        reduced_initial.push_back(std::move(reduced));
    }
    auto eig = detail::solve_lowest_sparse_modes(k, m, requested_modes, max_iterations,
                                                 tolerance, use_rcm_ordering, reduced_initial);
    auto count_positive = [](const detail::SparseEigenResult& candidate) {
        if (candidate.eigenvalues.empty()) return 0;
        const double scale = std::abs(candidate.eigenvalues.back());
        const double zero_tol = std::max(scale * 1e-12, std::numeric_limits<double>::min());
        int count = 0;
        for (double value : candidate.eigenvalues)
            if (value > zero_tol && std::isfinite(value)) ++count;
        return count;
    };
    if (eig.converged && count_positive(eig) < requested_modes) {
        const int retry_modes = requested_modes + std::max(2, requested_modes / 4);
        if (transform.reduced_dof_count > retry_modes) {
            auto retry = detail::solve_lowest_sparse_modes(
                k, m, retry_modes, max_iterations, tolerance, use_rcm_ordering, reduced_initial);
            if (retry.converged && count_positive(retry) >= requested_modes)
                eig = std::move(retry);
        }
    }
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
        for (int i = 0; i < n; ++i) {
            const int ri = transform.reduced_index[static_cast<std::size_t>(i)];
            if (ri >= 0) full[static_cast<std::size_t>(i)] =
                transform.scale[static_cast<std::size_t>(i)] * eig.vectors[mode][static_cast<std::size_t>(ri)];
        }
        out.eigenvalues.push_back(lambda);
        out.frequencies_hz.push_back(std::sqrt(lambda) / two_pi);
        out.mode_shapes.push_back(std::move(full));
        if (out.frequencies_hz.size() == static_cast<std::size_t>(requested_modes)) break;
    }
    if (out.frequencies_hz.size() < static_cast<std::size_t>(requested_modes)) out.converged = false;
    return out;
}

} // namespace BellOptimizer::FEM
