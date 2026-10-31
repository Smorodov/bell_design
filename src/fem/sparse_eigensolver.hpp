#pragma once

#include "sparse_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace BellOptimizer::detail {

struct SparseEigenResult {
    std::vector<double> eigenvalues;
    // Eigenvectors are stored column-wise: vectors[mode][dof].
    std::vector<std::vector<double>> vectors;
    int iterations = 0;
    bool converged = false;
    // Normwise backward error for K*x=lambda*M*x.
    double final_relative_residual = std::numeric_limits<double>::infinity();
    double final_eigenvalue_change = std::numeric_limits<double>::infinity();
    std::size_t ic0_nonzeros = 0;
    std::size_t pcg_iterations = 0;
    double max_mass_orthogonality_error = std::numeric_limits<double>::infinity();
    double max_eigenpair_relative_residual = std::numeric_limits<double>::infinity();
    double regularization_shift = std::numeric_limits<double>::quiet_NaN();
};

namespace sparse_eigen_detail {
inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
inline std::vector<double> add_scaled(const std::vector<double>& a,
                                      const std::vector<double>& b, double scale) {
    std::vector<double> out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) out[i] = a[i] + scale * b[i];
    return out;
}
inline bool finite_vector(const std::vector<double>& v) {
    for (double value : v) if (!std::isfinite(value)) return false;
    return true;
}

inline bool valid_csr(const SparseMatrix& a) {
    const std::size_t n = static_cast<std::size_t>(a.n);
    if (a.n <= 0 || a.row_offsets.size() != n + 1u ||
        a.row_offsets.front() != 0u || a.row_offsets.back() != a.columns.size() ||
        a.columns.size() != a.values.size()) return false;
    for (int row = 0; row < a.n; ++row) {
        const std::size_t begin = a.row_offsets[static_cast<std::size_t>(row)];
        const std::size_t end = a.row_offsets[static_cast<std::size_t>(row) + 1u];
        if (begin > end) return false;
        int previous_column = -1;
        for (std::size_t p = begin; p < end; ++p) {
            const int column = a.columns[p];
            if (column < 0 || column >= a.n || column <= previous_column ||
                !std::isfinite(a.values[p])) return false;
            previous_column = column;
        }
    }
    return true;
}

inline bool valid_symmetric_matrix(const SparseMatrix& a, bool require_positive_diagonal) {
    if (!valid_csr(a)) return false;
    constexpr double symmetry_tolerance = 1e-10;
    for (int row = 0; row < a.n; ++row) {
        const double diagonal = a.get(row, row);
        if (!std::isfinite(diagonal) ||
            (require_positive_diagonal && !(diagonal > 0.0))) return false;
        for (std::size_t p = a.row_offsets[static_cast<std::size_t>(row)];
             p < a.row_offsets[static_cast<std::size_t>(row) + 1u]; ++p) {
            const int column = a.columns[p];
            if (column == row) continue;
            const double aij = a.values[p];
            const double aji = a.get(column, row);
            const double scale = std::max({std::abs(aij), std::abs(aji),
                                           std::numeric_limits<double>::min()});
            // Check both triangular directions.  Iterating only over entries
            // above the diagonal would miss a one-sided lower-triangular
            // coefficient because the absent upper entry has no row to visit.
            if (std::abs(aij - aji) > symmetry_tolerance * scale) return false;
        }
    }
    return true;
}

struct IC0Preconditioner {
    int n = 0;
    std::vector<std::vector<std::pair<int, double>>> lower;
    // For the L^T solve, upper[column] contains (row, L(row,column)).
    std::vector<std::vector<std::pair<int, double>>> upper;
    std::vector<double> diag;
    bool valid = false;

    std::size_t factor_nonzeros() const {
        std::size_t count = diag.size();
        for (const auto& row : lower) count += row.size();
        return count;
    }

    bool build(const SparseMatrix& a, double diagonal_shift) {
        valid = false;
        n = a.n;
        lower.assign(static_cast<std::size_t>(n), {});
        upper.assign(static_cast<std::size_t>(n), {});
        diag.assign(static_cast<std::size_t>(n), 0.0);
        // IC(0) keeps exactly the lower-triangular sparsity pattern of A.
        // The diagonal shift is applied only to the factorization system;
        // it is not part of the physical eigenproblem.
        for (int i = 0; i < n; ++i) {
            for (std::size_t p = a.row_offsets[static_cast<std::size_t>(i)];
                 p < a.row_offsets[static_cast<std::size_t>(i) + 1]; ++p) {
                const int j = a.columns[p];
                const double value = a.values[p];
                if (j < i) lower[static_cast<std::size_t>(i)].push_back({j, value});
                else if (j == i) diag[static_cast<std::size_t>(i)] = value + diagonal_shift;
            }
        }
        for (int i = 0; i < n; ++i) {
            auto& row = lower[static_cast<std::size_t>(i)];
            std::sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
            for (auto& entry : row) {
                const int j = entry.first;
                double sum = entry.second;
                const auto& ri = row;
                const auto& rj = lower[static_cast<std::size_t>(j)];
                std::size_t pi = 0, pj = 0;
                while (pi < ri.size() && ri[pi].first < j && pj < rj.size()) {
                    if (ri[pi].first == rj[pj].first) {
                        sum -= ri[pi].second * rj[pj].second;
                        ++pi; ++pj;
                    } else if (ri[pi].first < rj[pj].first) ++pi;
                    else ++pj;
                }
                const double dj = diag[static_cast<std::size_t>(j)];
                if (!(dj > 0.0) || !std::isfinite(dj)) { valid = false; return false; }
                entry.second = sum / dj;
                upper[static_cast<std::size_t>(j)].push_back({i, entry.second});
            }
            double d = diag[static_cast<std::size_t>(i)];
            for (const auto& entry : row) d -= entry.second * entry.second;
            if (!(d > 0.0) || !std::isfinite(d)) { valid = false; return false; }
            diag[static_cast<std::size_t>(i)] = std::sqrt(d);
        }
        valid = true;
        return true;
    }

    bool solve(const std::vector<double>& r, std::vector<double>& z) const {
        if (!valid || static_cast<int>(r.size()) != n) return false;
        z = r;
        // L*y=r.
        for (int i = 0; i < n; ++i) {
            double sum = z[static_cast<std::size_t>(i)];
            for (const auto& e : lower[static_cast<std::size_t>(i)])
                sum -= e.second * z[static_cast<std::size_t>(e.first)];
            z[static_cast<std::size_t>(i)] = sum / diag[static_cast<std::size_t>(i)];
        }
        // L^T*z=y. Process each column directly through the transpose
        // adjacency, keeping the solve O(nnz(L)) instead of O(n^2).
        for (int i = n - 1; i >= 0; --i) {
            double sum = z[static_cast<std::size_t>(i)];
            for (const auto& e : upper[static_cast<std::size_t>(i)])
                sum -= e.second * z[static_cast<std::size_t>(e.first)];
            z[static_cast<std::size_t>(i)] = sum / diag[static_cast<std::size_t>(i)];
        }
        return true;
    }
};

inline bool pcg(const SparseMatrix& a, const std::vector<double>& b,
                std::vector<double>& x, const IC0Preconditioner& preconditioner,
                double tolerance, int max_iterations, std::size_t& iterations_used) {
    const std::size_t n = b.size();
    iterations_used = 0;
    if (n != static_cast<std::size_t>(a.n) || max_iterations <= 0 ||
        !(tolerance > 0.0) || !std::isfinite(tolerance) || !finite_vector(b)) return false;
    x.assign(n, 0.0);
    std::vector<double> r(b);
    std::vector<double> z(n, 0.0);
    std::vector<double> p(n, 0.0);
    if (!preconditioner.solve(r, z)) return false;
    p = z;
    double rz = dot(r, z);
    if (!(rz > 0.0) || !std::isfinite(rz)) return false;
    const double bnorm = std::sqrt(std::max(dot(b, b), std::numeric_limits<double>::min()));
    if (std::sqrt(dot(r, r)) <= tolerance * bnorm) return true;
    for (int it = 0; it < max_iterations; ++it) {
        iterations_used = static_cast<std::size_t>(it + 1);
        const auto ap = a.multiply(p);
        const double denominator = dot(p, ap);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) return false;
        const double alpha = rz / denominator;
        if (!std::isfinite(alpha)) return false;
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        const double residual_norm = std::sqrt(std::max(dot(r, r), 0.0));
        if (!std::isfinite(residual_norm)) return false;
        if (residual_norm <= tolerance * bnorm) return true;
        if (!preconditioner.solve(r, z)) return false;
        const double next_rz = dot(r, z);
        if (!(next_rz > 0.0) || !std::isfinite(next_rz)) return false;
        const double beta = next_rz / rz;
        if (!(beta >= 0.0) || !std::isfinite(beta)) return false;
        for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz = next_rz;
    }
    return false;
}

inline bool jacobi(std::vector<double> a, int n, std::vector<double>& eval,
                   std::vector<double>& evec) {
    evec.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) evec[static_cast<std::size_t>(i) * n + i] = 1.0;
    const int limit = std::max(100, 80 * n * n);
    bool converged = false;
    double frobenius2 = 0.0;
    for (double value : a) {
        if (!std::isfinite(value)) return false;
        frobenius2 += value * value;
    }
    if (!(frobenius2 >= 0.0) || !std::isfinite(frobenius2)) return false;
    const double matrix_scale = std::max(std::sqrt(frobenius2),
                                         std::numeric_limits<double>::denorm_min());
    const double off_diagonal_tolerance = 1e-14 * matrix_scale;
    for (int it = 0; it < limit; ++it) {
        int p = 0, q = 0; double largest = 0.0;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
            const double v = std::abs(a[static_cast<std::size_t>(i) * n + j]);
            if (v > largest) { largest = v; p = i; q = j; }
        }
        if (largest <= off_diagonal_tolerance) { converged = true; break; }
        const auto at = [n](int i, int j) { return static_cast<std::size_t>(i) * n + j; };
        const double app = a[at(p,p)], aqq = a[at(q,q)], apq = a[at(p,q)];
        const double angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
        const double c = std::cos(angle), s = std::sin(angle);
        for (int k = 0; k < n; ++k) {
            const double x = a[at(p,k)], y = a[at(q,k)];
            a[at(p,k)] = c*x - s*y; a[at(q,k)] = s*x + c*y;
        }
        for (int k = 0; k < n; ++k) {
            const double x = a[at(k,p)], y = a[at(k,q)];
            a[at(k,p)] = c*x - s*y; a[at(k,q)] = s*x + c*y;
            const double vx = evec[at(k,p)], vy = evec[at(k,q)];
            evec[at(k,p)] = c*vx - s*vy; evec[at(k,q)] = s*vx + c*vy;
        }
    }
    if (!converged) return false;
    eval.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) eval[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i) * n + i];
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
        if (eval[static_cast<std::size_t>(j)] < eval[static_cast<std::size_t>(i)]) {
            std::swap(eval[static_cast<std::size_t>(i)], eval[static_cast<std::size_t>(j)]);
            for (int k = 0; k < n; ++k) std::swap(evec[static_cast<std::size_t>(k)*n+i], evec[static_cast<std::size_t>(k)*n+j]);
        }
    return true;
}
} // namespace sparse_eigen_detail

// Dependency-free block inverse/subspace iteration for the lowest generalized
// symmetric modes K*x=lambda*M*x. Uses a shifted sparse CG solve and a small
// dense Rayleigh-Ritz problem. Intended as a first-mode solver; CG failure or
// loss of subspace rank is reported instead of silently returning modes.
inline SparseEigenResult solve_lowest_sparse_modes(const SparseMatrix& k, const SparseMatrix& m,
                                                    int requested_modes, int max_iterations = 120,
                                                    double tolerance = 1e-8, bool use_rcm_ordering = false,
                                                    const std::vector<std::vector<double>>& initial_vectors = {}) {
    using namespace sparse_eigen_detail;
    SparseEigenResult result;
    if (k.n <= 0 || k.n != m.n || requested_modes <= 0 || max_iterations <= 0 ||
        !(tolerance > 0.0) || !std::isfinite(tolerance) ||
        !valid_symmetric_matrix(k, false) || !valid_symmetric_matrix(m, true)) return result;
    const int n = k.n;
    if (n > (std::numeric_limits<int>::max() / 8)) return result;
    if (use_rcm_ordering) {
        // One common permutation must be applied to both matrices of the
        // generalized eigenproblem. Solve entirely in the permuted basis and
        // map the eigenvectors back afterwards. The physical matrices are
        // unchanged up to relabeling; no coefficients are dropped.
        const auto graph = structural_union(k, m);
        const auto order = reverse_cuthill_mckee_order(graph);
        const auto kp = symmetric_permute(k, order);
        const auto mp = symmetric_permute(m, order);

        // symmetric_permute() uses order[new] = old. Initial vectors are
        // physical vectors in the original DOF numbering, so they must be
        // transformed into the same permuted basis before the recursive solve.
        // Passing them unchanged seeds the solver with values attached to the
        // wrong DOFs and breaks continuation when RCM ordering is enabled.
        std::vector<std::vector<double>> permuted_initial;
        permuted_initial.reserve(initial_vectors.size());
        for (const auto& vector : initial_vectors) {
            if (vector.size() != static_cast<std::size_t>(n)) continue;
            std::vector<double> reordered(static_cast<std::size_t>(n), 0.0);
            for (int new_index = 0; new_index < n; ++new_index)
                reordered[static_cast<std::size_t>(new_index)] =
                    vector[static_cast<std::size_t>(order[static_cast<std::size_t>(new_index)])];
            permuted_initial.push_back(std::move(reordered));
        }

        auto permuted = solve_lowest_sparse_modes(kp, mp, requested_modes,
                                                  max_iterations, tolerance, false, permuted_initial);
        for (auto& vector : permuted.vectors) {
            std::vector<double> original(static_cast<std::size_t>(n), 0.0);
            for (int new_index = 0; new_index < n; ++new_index)
                original[static_cast<std::size_t>(order[static_cast<std::size_t>(new_index)])] =
                    vector[static_cast<std::size_t>(new_index)];
            vector.swap(original);
        }
        return permuted;
    }
    const int block = std::min(n, requested_modes + std::max(3, requested_modes / 2));
    double max_k_diag = 0.0, max_m_diag = 0.0;
    for (int i = 0; i < n; ++i) {
        max_k_diag = std::max(max_k_diag, std::abs(k.get(i, i)));
        max_m_diag = std::max(max_m_diag, std::abs(m.get(i, i)));
    }
    if (!(max_k_diag > 0.0) || !(max_m_diag > 0.0)) return result;
    // Normalize K and M independently before inverse iteration.  The
    // generalized problem is invariant under K -> K/k_scale and
    // M -> M/m_scale; the physical eigenvalue is restored by
    // lambda = lambda_normalized * (k_scale/m_scale).  This prevents very
    // large or very small absolute matrix coefficients from degrading PCG
    // and the small dense Ritz problem.
    const double k_scale = max_k_diag;
    const double m_scale = max_m_diag;
    const double spectral_scale = k_scale / m_scale;
    if (!(spectral_scale > 0.0) || !std::isfinite(spectral_scale)) return result;

    std::vector<SparseTriplet> normalized_k_entries;
    std::vector<SparseTriplet> normalized_m_entries;
    normalized_k_entries.reserve(k.values.size());
    normalized_m_entries.reserve(m.values.size());
    for (int i = 0; i < n; ++i) {
        for (std::size_t p = k.row_offsets[static_cast<std::size_t>(i)];
             p < k.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            normalized_k_entries.push_back({i, k.columns[p], k.values[p] / k_scale});
        for (std::size_t p = m.row_offsets[static_cast<std::size_t>(i)];
             p < m.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            normalized_m_entries.push_back({i, m.columns[p], m.values[p] / m_scale});
    }
    const SparseMatrix normalized_k = SparseMatrix::from_triplets(n, std::move(normalized_k_entries));
    const SparseMatrix normalized_m = SparseMatrix::from_triplets(n, std::move(normalized_m_entries));

    double k_frobenius2 = 0.0, m_frobenius2 = 0.0;
    for (double value : normalized_k.values) k_frobenius2 += value * value;
    for (double value : normalized_m.values) m_frobenius2 += value * value;
    const double k_norm = std::sqrt(k_frobenius2);
    const double m_norm = std::sqrt(m_frobenius2);
    if (!(k_norm > 0.0) || !(m_norm > 0.0) || !std::isfinite(k_norm) || !std::isfinite(m_norm)) return result;

    // A positive normalized shift regularizes a semidefinite/near-null
    // stiffness system.  It affects only the linear solves, never the Ritz
    // values reconstructed from the physical K/M pair.
    const double shift = 1e-12;
    IC0Preconditioner preconditioner;
    double factor_shift = shift;
    bool factor_ok = false;
    for (int attempt = 0; attempt < 12 && !factor_ok; ++attempt) {
        std::vector<SparseTriplet> regularized_entries;
        regularized_entries.reserve(normalized_k.values.size() + normalized_m.values.size());
        for (int i = 0; i < n; ++i) {
            for (std::size_t p = normalized_k.row_offsets[static_cast<std::size_t>(i)];
                 p < normalized_k.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
                regularized_entries.push_back({i, normalized_k.columns[p], normalized_k.values[p]});
            for (std::size_t p = normalized_m.row_offsets[static_cast<std::size_t>(i)];
                 p < normalized_m.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
                regularized_entries.push_back({i, normalized_m.columns[p], factor_shift * normalized_m.values[p]});
        }
        const SparseMatrix regularized = SparseMatrix::from_triplets(n, std::move(regularized_entries));
        factor_ok = preconditioner.build(regularized, 0.0);
        if (!factor_ok) factor_shift *= 10.0;
    }
    if (!factor_ok) return result;
    result.ic0_nonzeros = preconditioner.factor_nonzeros();
    result.regularization_shift = factor_shift;
    std::vector<SparseTriplet> final_entries;
    final_entries.reserve(normalized_k.values.size() + normalized_m.values.size());
    for (int i = 0; i < n; ++i) {
        for (std::size_t p = normalized_k.row_offsets[static_cast<std::size_t>(i)];
             p < normalized_k.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            final_entries.push_back({i, normalized_k.columns[p], normalized_k.values[p]});
        for (std::size_t p = normalized_m.row_offsets[static_cast<std::size_t>(i)];
             p < normalized_m.row_offsets[static_cast<std::size_t>(i) + 1u]; ++p)
            final_entries.push_back({i, normalized_m.columns[p], factor_shift * normalized_m.values[p]});
    }
    const SparseMatrix solve_matrix = SparseMatrix::from_triplets(n, std::move(final_entries));
    std::vector<std::vector<double>> q(static_cast<std::size_t>(block), std::vector<double>(static_cast<std::size_t>(n)));
    if (!initial_vectors.empty()) {
        const std::size_t seeded = std::min(initial_vectors.size(), q.size());
        for (std::size_t j = 0; j < seeded; ++j) {
            if (initial_vectors[j].size() != static_cast<std::size_t>(n)) return result;
            q[j] = initial_vectors[j];
        }
        for (std::size_t j = seeded; j < q.size(); ++j)
            for (int i = 0; i < n; ++i)
                q[j][static_cast<std::size_t>(i)] =
                    std::sin((i + 1.0) * (j + 1.0) * 0.731) + std::cos((i + 2.0) * (j + 1.0) * 0.193);
    } else {
        for (int j = 0; j < block; ++j) for (int i = 0; i < n; ++i)
            q[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] =
                std::sin((i + 1.0) * (j + 1.0) * 0.731) + std::cos((i + 2.0) * (j + 1.0) * 0.193);
    }
    auto m_orthonormalize = [&](std::vector<std::vector<double>>& basis) {
        int kept = 0;
        for (std::size_t source = 0; source < basis.size(); ++source) {
            auto& v = basis[source];
            if (v.size() != static_cast<std::size_t>(n) || !finite_vector(v)) continue;

            // Two-pass modified Gram-Schmidt is used here deliberately.  The
            // first pass removes the dominant projection; the second pass
            // restores orthogonality lost to cancellation, which matters for
            // clustered generalized eigenvalues and continuation vectors.
            for (int pass = 0; pass < 2; ++pass) {
                auto mv = normalized_m.multiply(v);
                for (int j = 0; j < kept; ++j) {
                    const auto& basis_j = basis[static_cast<std::size_t>(j)];
                    const double c = dot(basis_j, mv);
                    if (!std::isfinite(c)) {
                        v.clear();
                        break;
                    }
                    const auto mbasis_j = normalized_m.multiply(basis_j);
                    for (int i = 0; i < n; ++i) {
                        v[static_cast<std::size_t>(i)] -= c * basis_j[static_cast<std::size_t>(i)];
                        mv[static_cast<std::size_t>(i)] -= c * mbasis_j[static_cast<std::size_t>(i)];
                    }
                }
                if (v.empty()) break;
            }
            if (v.empty() || !finite_vector(v)) continue;

            const auto mv2 = normalized_m.multiply(v);
            const double norm2 = dot(v, mv2);
            if (!(norm2 > 1e-24) || !std::isfinite(norm2)) continue;
            const double inv = 1.0 / std::sqrt(norm2);
            if (!std::isfinite(inv)) continue;
            for (double& x : v) x *= inv;
            if (!finite_vector(v)) continue;
            if (static_cast<std::size_t>(kept) != source)
                basis[static_cast<std::size_t>(kept)] = std::move(v);
            ++kept;
        }
        basis.resize(static_cast<std::size_t>(kept));
        return kept;
    };
    if (m_orthonormalize(q) < std::min(block, requested_modes)) return result;
    std::vector<double> previous(static_cast<std::size_t>(requested_modes), std::numeric_limits<double>::infinity());
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<std::vector<double>> next;
        next.reserve(q.size());
        for (const auto& v : q) {
            const auto rhs = normalized_m.multiply(v);
            std::vector<double> x;
            std::size_t pcg_used = 0;
            const double pcg_tolerance = std::max(1e-14, std::min(1e-10, 0.1 * tolerance));
            if (!pcg(solve_matrix, rhs, x, preconditioner, pcg_tolerance, std::max(1000, 32 * n), pcg_used)) return result;
            result.pcg_iterations += pcg_used;
            next.push_back(std::move(x));
        }
        if (m_orthonormalize(next) < requested_modes) return result;
        q = std::move(next);
        const int b = static_cast<int>(q.size());
        std::vector<double> projected_k(static_cast<std::size_t>(b)*b, 0.0);
        for (int j = 0; j < b; ++j) {
            const auto kj = normalized_k.multiply(q[static_cast<std::size_t>(j)]);
            for (int i = 0; i < b; ++i)
                projected_k[static_cast<std::size_t>(i)*b+j] = dot(q[static_cast<std::size_t>(i)], kj);
        }
        // M-orthonormal basis makes projected M identity to roundoff.
        // The dense Jacobi stage must preserve that metric after rotating the
        // Ritz basis; verify the returned Ritz vectors explicitly instead of
        // assuming the projected algebra was exact.
        std::vector<double> eval, evec;
        if (!jacobi(projected_k, b, eval, evec)) return result;
        result.eigenvalues.clear(); result.vectors.clear();
        double max_change = 0.0;
        double projected_spectral_scale = 0.0;
        for (double value : eval) projected_spectral_scale = std::max(projected_spectral_scale, std::abs(value));
        const double convergence_floor = std::max(projected_spectral_scale * 1e-8, 1e-30);
        for (int mode = 0; mode < std::min(requested_modes, b); ++mode) {
            std::vector<double> v(static_cast<std::size_t>(n), 0.0);
            for (int j = 0; j < b; ++j) {
                const double c = evec[static_cast<std::size_t>(j)*b+mode];
                for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] += c*q[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
            }
            const double physical_eigenvalue = eval[static_cast<std::size_t>(mode)] * spectral_scale;
            result.eigenvalues.push_back(physical_eigenvalue);
            result.vectors.push_back(std::move(v));
            const double change = std::isfinite(previous[static_cast<std::size_t>(mode)])
                ? std::abs(eval[static_cast<std::size_t>(mode)] - previous[static_cast<std::size_t>(mode)]) /
                  std::max(std::abs(eval[static_cast<std::size_t>(mode)]), convergence_floor)
                : std::numeric_limits<double>::infinity();
            max_change = std::max(max_change, change);
            previous[static_cast<std::size_t>(mode)] = eval[static_cast<std::size_t>(mode)];
        }
        result.iterations = iteration + 1;
        // Eigenvalue stagnation alone is not sufficient for clustered modes:
        // their values may stabilize while individual Ritz vectors still mix.
        // Require a scale-aware generalized eigenpair residual as well.
        double max_relative_residual = 0.0;
        for (std::size_t mode = 0; mode < result.vectors.size(); ++mode) {
            const double normalized_eigenvalue = eval[mode];
            const auto kv = normalized_k.multiply(result.vectors[mode]);
            const auto mv = normalized_m.multiply(result.vectors[mode]);
            double residual2 = 0.0, vector_norm2 = 0.0;
            for (int i = 0; i < n; ++i) {
                const double lambda_mv = normalized_eigenvalue * mv[static_cast<std::size_t>(i)];
                const double residual = kv[static_cast<std::size_t>(i)] - lambda_mv;
                residual2 += residual * residual;
                vector_norm2 += result.vectors[mode][static_cast<std::size_t>(i)] *
                                result.vectors[mode][static_cast<std::size_t>(i)];
            }
            // Normwise backward error for K*x=lambda*M*x. Unlike the former
            // ||r||/(||Kx||+|lambda| ||Mx||) scale, this remains meaningful when
            // both K*x and lambda*M*x are close to zero.
            const double matrix_scale = k_norm + std::abs(normalized_eigenvalue) * m_norm;
            const double scale = std::max(matrix_scale * std::sqrt(vector_norm2),
                                          std::numeric_limits<double>::min());
            max_relative_residual = std::max(max_relative_residual, std::sqrt(residual2) / scale);
        }

        double max_mass_orthogonality_error = 0.0;
        for (std::size_t i = 0; i < result.vectors.size(); ++i) {
            const auto mi = normalized_m.multiply(result.vectors[i]);
            for (std::size_t j = 0; j < result.vectors.size(); ++j) {
                const double gij = dot(result.vectors[j], mi);
                const double target = (i == j) ? 1.0 : 0.0;
                if (!std::isfinite(gij)) return result;
                max_mass_orthogonality_error =
                    std::max(max_mass_orthogonality_error, std::abs(gij - target));
            }
        }
        result.max_eigenpair_relative_residual = max_relative_residual;
        result.max_mass_orthogonality_error = max_mass_orthogonality_error;
        result.final_relative_residual = max_relative_residual;
        result.final_eigenvalue_change = max_change;
        const double residual_tolerance = std::max(10.0 * tolerance, 1e-10);
        const double orthogonality_tolerance = 1e-8;
        if (max_change < tolerance && max_relative_residual < residual_tolerance &&
            max_mass_orthogonality_error < orthogonality_tolerance) {
            result.converged = true;
            return result;
        }
    }
    return result;
}

} // namespace BellOptimizer::detail
