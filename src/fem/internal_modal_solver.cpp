#include "internal_modal_solver.hpp"
#include "sparse_matrix.hpp"
#include "sparse_eigensolver.hpp"
#include "bell_types.h"
#include "math_constants.hpp"
#include "modal_influence.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

namespace BellOptimizer {
namespace {

constexpr double PI = bell::kPi;

struct Matrix {
    int n = 0;
    std::vector<double> a;
    explicit Matrix(int n_=0) : n(n_), a(static_cast<size_t>(n_) * static_cast<size_t>(n_), 0.0) {}
    double& operator()(int i, int j) { return a[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)]; }
    double operator()(int i, int j) const { return a[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)]; }
};

using Triplet = detail::SparseTriplet;
using SparseMatrix = detail::SparseMatrix;

Matrix sparse_to_dense(const SparseMatrix& sparse) {
    Matrix dense(sparse.n);
    for (int row = 0; row < sparse.n; ++row)
        for (size_t k = sparse.row_offsets[static_cast<size_t>(row)];
             k < sparse.row_offsets[static_cast<size_t>(row) + 1u]; ++k)
            dense(row, sparse.columns[k]) = sparse.values[k];
    return dense;
}

MatrixAuditStats audit_matrix(const SparseMatrix& A) {
    MatrixAuditStats stats;
    stats.dimension = static_cast<size_t>(A.n);
    const size_t n = static_cast<size_t>(A.n);
    stats.dense_bytes = n * n * sizeof(double);
    stats.exact_nonzeros = A.values.size();
    stats.density = n == 0 ? 0.0 : static_cast<double>(stats.exact_nonzeros) /
        static_cast<double>(n * n);
    stats.estimated_csr_bytes = A.storage_bytes();

    double norm_sq = 0.0;
    double asymmetry_sq = 0.0;
    double min_positive_diag = std::numeric_limits<double>::infinity();
    double max_positive_diag = 0.0;
    size_t nonpositive_diag = 0;
    for (int i = 0; i < A.n; ++i) {
        const double diagonal = A.get(i, i);
        if (diagonal > 0.0 && std::isfinite(diagonal)) {
            min_positive_diag = std::min(min_positive_diag, diagonal);
            max_positive_diag = std::max(max_positive_diag, diagonal);
        } else {
            ++nonpositive_diag;
        }
    }
    if (std::isfinite(min_positive_diag) && max_positive_diag > 0.0) {
        stats.positive_diagonal_min = min_positive_diag;
        stats.positive_diagonal_max = max_positive_diag;
        stats.positive_diagonal_dynamic_range = max_positive_diag / min_positive_diag;
    }
    stats.nonpositive_diagonal_count = nonpositive_diag;
    std::size_t total_row_nnz = 0;
    std::size_t max_bandwidth = 0;
    std::size_t matrix_profile = 0;
    for (int row = 0; row < A.n; ++row) {
        const std::size_t begin = A.row_offsets[static_cast<std::size_t>(row)];
        const std::size_t end = A.row_offsets[static_cast<std::size_t>(row) + 1u];
        total_row_nnz += end - begin;
        std::size_t first_column = static_cast<std::size_t>(row);
        bool has_offdiagonal = false;
        for (std::size_t p = begin; p < end; ++p) {
            const std::size_t column = static_cast<std::size_t>(A.columns[p]);
            if (column != static_cast<std::size_t>(row)) {
                has_offdiagonal = true;
                max_bandwidth = std::max(max_bandwidth,
                    column > static_cast<std::size_t>(row)
                        ? column - static_cast<std::size_t>(row)
                        : static_cast<std::size_t>(row) - column);
                first_column = std::min(first_column, column);
            }
        }
        if (has_offdiagonal)
            matrix_profile += static_cast<std::size_t>(row) - first_column;
    }
    stats.mean_row_nnz = n == 0 ? 0.0 : static_cast<double>(total_row_nnz) / static_cast<double>(n);
    stats.max_bandwidth = max_bandwidth;
    stats.profile = matrix_profile;

    // RCM is diagnostic-only at this stage. The production ordering must be
    // common to K and M, so the actual solve is not changed here. The common
    // FEM graph is installed below by the solve() audit after both matrices
    // are available; these per-matrix fields are populated there.

    std::vector<double> row_scale(n, 0.0);
    for (int row = 0; row < A.n; ++row) {
        for (size_t k = A.row_offsets[static_cast<size_t>(row)];
             k < A.row_offsets[static_cast<size_t>(row) + 1u]; ++k) {
            const double value = A.values[k];
            norm_sq += value * value;
            row_scale[static_cast<size_t>(row)] =
                std::max(row_scale[static_cast<size_t>(row)], std::abs(value));
            const int column = A.columns[k];
            if (column > row) {
                const double difference = value - A.get(column, row);
                asymmetry_sq += 2.0 * difference * difference;
                const double pair_magnitude = std::max(std::abs(value), std::abs(A.get(column, row)));
                const double local_scale = std::max(row_scale[static_cast<size_t>(row)],
                                                    row_scale[static_cast<size_t>(column)]);
                // Row scales are completed in a second pass below. Count thresholds there.
                (void)pair_magnitude;
                (void)local_scale;
            }
        }
    }
    stats.relative_asymmetry_frobenius = std::sqrt(asymmetry_sq) /
        std::max(std::sqrt(norm_sq), std::numeric_limits<double>::min());
    // Structural pattern audit: FEM assembly should produce a symmetric graph.
    // Count only stored off-diagonal pairs; this is independent of coefficient size.
    for (int row = 0; row < A.n; ++row) {
        for (std::size_t k = A.row_offsets[static_cast<std::size_t>(row)];
             k < A.row_offsets[static_cast<std::size_t>(row) + 1u]; ++k) {
            const int column = A.columns[k];
            if (column <= row) continue;
            if (A.get(column, row) != 0.0) ++stats.symmetric_offdiagonal_pairs;
            else ++stats.structural_asymmetry_pairs;
        }
    }
    // Threshold diagnostics only. Scan stored upper-triangle pairs after all row
    // scales have been computed; structural zeros are not candidate removals.
    for (int row = 0; row < A.n; ++row) {
        for (size_t k = A.row_offsets[static_cast<size_t>(row)];
             k < A.row_offsets[static_cast<size_t>(row) + 1u]; ++k) {
            const int column = A.columns[k];
            if (column <= row) continue;
            const double pair_magnitude = std::max(std::abs(A.values[k]), std::abs(A.get(column, row)));
            const double local_scale = std::max(row_scale[static_cast<size_t>(row)],
                                                row_scale[static_cast<size_t>(column)]);
            if (pair_magnitude == 0.0 || local_scale == 0.0) continue;
            const double ratio = pair_magnitude / local_scale;
            if (ratio <= 1e-12) ++stats.candidate_pairs_below_1e12;
            if (ratio <= 1e-10) ++stats.candidate_pairs_below_1e10;
            if (ratio <= 1e-8) ++stats.candidate_pairs_below_1e8;
        }
    }
    return stats;
}

bool cholesky(const Matrix& M, Matrix& L) {
    const int n = M.n;
    L = Matrix(n);
    for (int i=0;i<n;++i) {
        for (int j=0;j<=i;++j) {
            double s = M(i,j);
            for (int k=0;k<j;++k) s -= L(i,k)*L(j,k);
            if (i == j) {
                if (!(s > 1e-18) || !std::isfinite(s)) return false;
                L(i,j) = std::sqrt(s);
            } else {
                L(i,j) = s / L(j,j);
            }
        }
    }
    return true;
}

// Symmetric Jacobi eigensolver. Small dense matrices are deliberate here:
// this is the validation implementation, not the final large-scale solver.
bool jacobi(Matrix A, std::vector<double>& eval, Matrix& evec) {
    const int n = A.n;
    evec = Matrix(n);
    for (int i=0;i<n;++i) evec(i,i)=1.0;
    const int max_iter = std::max(100, 80*n*n);
    bool converged = false;
    for (int it=0; it<max_iter; ++it) {
        int p=0,q=1;
        double mx=0.0;
        for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
            const double x=std::abs(A(i,j));
            if (x>mx) { mx=x; p=i; q=j; }
        }
        if (mx < 1e-11) { converged = true; break; }
        const double app=A(p,p), aqq=A(q,q), apq=A(p,q);
        const double phi=0.5*std::atan2(2.0*apq, aqq-app);
        const double c=std::cos(phi), s=std::sin(phi);
        for (int k=0;k<n;++k) {
            const double aik=A(p,k), aqk=A(q,k);
            A(p,k)=c*aik-s*aqk;
            A(q,k)=s*aik+c*aqk;
        }
        for (int k=0;k<n;++k) {
            const double akp=A(k,p), akq=A(k,q);
            A(k,p)=c*akp-s*akq;
            A(k,q)=s*akp+c*akq;
        }
        for (int k=0;k<n;++k) {
            const double vip=evec(k,p), viq=evec(k,q);
            evec(k,p)=c*vip-s*viq;
            evec(k,q)=s*vip+c*viq;
        }
    }
    // Never return diagonal entries as eigenvalues if the iteration budget was
    // exhausted. They would look plausible while eigenvectors may still be wrong.
    if (!converged) return false;
    eval.resize(n);
    for (int i=0;i<n;++i) eval[i]=A(i,i);
    for (int i=0;i<n;++i) for (int j=i+1;j<n;++j)
        if (eval[j] < eval[i]) {
            std::swap(eval[i],eval[j]);
            for (int k=0;k<n;++k) std::swap(evec(k,i),evec(k,j));
        }
    return true;
}

// Solve L X = B for each column, then A = L^-1 K L^-T.
void solve_upper_transpose(const Matrix& L, const std::vector<double>& b, std::vector<double>& x) {
    const int n = L.n;
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = b[static_cast<size_t>(i)];
        for (int k = i + 1; k < n; ++k) s -= L(k, i) * x[static_cast<size_t>(k)];
        x[static_cast<size_t>(i)] = s / L(i, i);
    }
}

Matrix generalized_to_standard(const Matrix& K, const Matrix& L) {
    const int n=K.n;
    Matrix X(n), A(n);
    for (int j=0;j<n;++j) {
        for (int i=0;i<n;++i) {
            double s=K(i,j);
            for (int k=0;k<i;++k) s-=L(i,k)*X(k,j);
            X(i,j)=s/L(i,i);
        }
    }
    // A = X * L^-T. Back-substitute by columns of L^T.
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) {
            double s=X(i,j);
            for (int k=0;k<j;++k) s-=A(i,k)*L(j,k);
            A(i,j)=s/L(j,j);
        }
    }
    // Symmetrize accumulated roundoff.
    for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
        const double v=0.5*(A(i,j)+A(j,i)); A(i,j)=A(j,i)=v;
    }
    return A;
}

bool element_matrices(const std::array<Node,4>& nd, double E, double nu, double rho,
                      double ke[8][8], double me[8][8]) {
    for(int i=0;i<8;++i) for(int j=0;j<8;++j) ke[i][j]=me[i][j]=0.0;
    const double lam=E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double mu=E/(2.0*(1.0+nu));
    double D[4][4] = {
        {lam+2*mu, lam,       lam,       0},
        {lam,       lam+2*mu, lam,       0},
        {lam,       lam,       lam+2*mu, 0},
        {0,         0,         0,         mu}
    };
    const double g=1.0/std::sqrt(3.0);
    for (double xi : {-g,g}) for (double eta : {-g,g}) {
        const double N[4] = {
            0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta),
            0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta)
        };
        const double dNdxi[4] = {-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)};
        const double dNdeta[4] = {-0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi)};
        double J00=0,J01=0,J10=0,J11=0;
        for(int a=0;a<4;++a){ J00+=dNdxi[a]*nd[a].r; J01+=dNdxi[a]*nd[a].z; J10+=dNdeta[a]*nd[a].r; J11+=dNdeta[a]*nd[a].z; }
        const double det=J00*J11-J01*J10;
        if (!(det>1e-14) || !std::isfinite(det)) return false;
        const double inv00= J11/det, inv01=-J01/det, inv10=-J10/det, inv11=J00/det;
        double dNdr[4], dNdz[4];
        for(int a=0;a<4;++a){ dNdr[a]=dNdxi[a]*inv00+dNdeta[a]*inv01; dNdz[a]=dNdxi[a]*inv10+dNdeta[a]*inv11; }
        double r=0; for(int a=0;a<4;++a) r+=N[a]*nd[a].r;
        // Gauss points must lie strictly inside the physical half-plane r > 0.
        // Clamping a non-positive radius here would hide an invalid/cross-axis
        // element and make the hoop strain inconsistent with the integration weight.
        if (!(r > 0.0) || !std::isfinite(r)) return false;
        const double r_safe=r;
        double B[4][8]{};
        for(int a=0;a<4;++a){
            const int u=2*a, v=u+1;
            B[0][u]=dNdr[a];
            B[1][v]=dNdz[a];
            B[2][u]=N[a]/r_safe;
            B[3][u]=dNdz[a]; B[3][v]=dNdr[a];
        }
        const double w=2.0*PI*r*det;
        for(int i=0;i<8;++i) for(int j=0;j<8;++j){
            double s=0; for(int a=0;a<4;++a) for(int b=0;b<4;++b) s+=B[a][i]*D[a][b]*B[b][j];
            ke[i][j]+=s*w;
        }
        for(int a=0;a<4;++a) for(int b=0;b<4;++b){
            const double mm=2.0*PI*r*det*rho*N[a]*N[b];
            me[2*a][2*b]+=mm; me[2*a+1][2*b+1]+=mm;
        }
    }
    return true;
}

} // namespace

bool InternalModalSolver::solve(const ModalProblem& problem, ModalData& result) {
    result.clear();
    const auto& mesh = problem.mesh;
    const auto& material = problem.material;
    if (mesh.nodes.empty() || mesh.elements.empty() ||
        mesh.nodes.size() > static_cast<size_t>(std::numeric_limits<int>::max() / 2) ||
        !std::isfinite(material.E) || !std::isfinite(material.rho) || !std::isfinite(material.nu) ||
        material.E <= 0.0 || material.rho <= 0.0 || material.nu <= -1.0 || material.nu >= 0.5) {
        if (!quiet_) std::cerr << "Internal FEM: invalid mesh/material or unsupported mesh size\n";
        return false;
    }

    const int node_count = static_cast<int>(mesh.nodes.size());
    const int ndof_full = node_count * 2;
    std::vector<int> id(mesh.nodes.size()+1,-1);
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const int node_id = mesh.nodes[i].id;
        const auto& node = mesh.nodes[i];
        if (node_id <= 0 || node_id > node_count || id[static_cast<size_t>(node_id)] != -1 ||
            !std::isfinite(node.r) || !std::isfinite(node.z) || node.r < -1e-9) {
            if (!quiet_) std::cerr << "Internal FEM: invalid, duplicate, or non-contiguous node IDs/coordinates\n";
            return false;
        }
        id[static_cast<size_t>(node_id)] = static_cast<int>(i);
    }

    // Apply the axis constraint before assembly. Constrained radial DOFs at r=0
    // never enter the global sparse matrices.
    std::vector<int> free_dof;
    free_dof.reserve(static_cast<size_t>(ndof_full));
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (std::abs(mesh.nodes[i].r) < 1e-9) free_dof.push_back(static_cast<int>(2u*i + 1u));
        else {
            free_dof.push_back(static_cast<int>(2u*i));
            free_dof.push_back(static_cast<int>(2u*i + 1u));
        }
    }
    const int n = static_cast<int>(free_dof.size());
    std::vector<int> full_to_reduced(static_cast<size_t>(ndof_full), -1);
    for (int reduced = 0; reduced < n; ++reduced)
        full_to_reduced[static_cast<size_t>(free_dof[static_cast<size_t>(reduced)])] = reduced;

    std::vector<Triplet> k_entries, m_entries;
    // Independent global accumulation path used only for small audit meshes.
    // Local element matrices are shared, but dense accumulation and CSR duplicate
    // reduction are separate; this checks DOF mapping and global assembly storage.
    constexpr int dense_reference_limit = 128;
    const bool check_dense_reference = n <= dense_reference_limit;
    std::vector<double> dense_reference_k, dense_reference_m;
    if (check_dense_reference) {
        const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
        dense_reference_k.assign(nn, 0.0);
        dense_reference_m.assign(nn, 0.0);
    }
    if (mesh.elements.size() > std::numeric_limits<size_t>::max() / 64u) {
        if (!quiet_) std::cerr << "Internal FEM: element count exceeds sparse assembly limits\n";
        return false;
    }
    const size_t triplet_reserve = mesh.elements.size() * 64u;
    k_entries.reserve(triplet_reserve);
    m_entries.reserve(triplet_reserve);
    for (const auto& e : mesh.elements) {
        if (e.n1 <= 0 || e.n2 <= 0 || e.n3 <= 0 || e.n4 <= 0 ||
            e.n1 > node_count || e.n2 > node_count || e.n3 > node_count || e.n4 > node_count ||
            id[static_cast<size_t>(e.n1)] < 0 || id[static_cast<size_t>(e.n2)] < 0 ||
            id[static_cast<size_t>(e.n3)] < 0 || id[static_cast<size_t>(e.n4)] < 0) {
            if (!quiet_) std::cerr << "Internal FEM: element references an invalid node ID\n";
            return false;
        }
        const std::array<Node,4> nd={mesh.nodes[static_cast<size_t>(id[static_cast<size_t>(e.n1)])],
                                     mesh.nodes[static_cast<size_t>(id[static_cast<size_t>(e.n2)])],
                                     mesh.nodes[static_cast<size_t>(id[static_cast<size_t>(e.n3)])],
                                     mesh.nodes[static_cast<size_t>(id[static_cast<size_t>(e.n4)])]};
        double ke[8][8], me[8][8];
        if (!element_matrices(nd,material.E,material.nu,material.rho,ke,me)) {
            if (!quiet_) std::cerr << "Internal FEM: invalid Q4 Jacobian\n";
            return false;
        }
        const int ids[4]={id[static_cast<size_t>(e.n1)],id[static_cast<size_t>(e.n2)],
                          id[static_cast<size_t>(e.n3)],id[static_cast<size_t>(e.n4)]};
        int element_dofs[8];
        for (int a = 0; a < 4; ++a) {
            element_dofs[2*a] = 2*ids[a];
            element_dofs[2*a+1] = 2*ids[a]+1;
        }
        for (int a = 0; a < 8; ++a) {
            const int row = full_to_reduced[static_cast<size_t>(element_dofs[a])];
            if (row < 0) continue;
            for (int b = 0; b < 8; ++b) {
                const int column = full_to_reduced[static_cast<size_t>(element_dofs[b])];
                if (column < 0) continue;
                const double kv = ke[a][b];
                const double mv = me[a][b];
                if (kv != 0.0) k_entries.push_back({row, column, kv});
                if (mv != 0.0) m_entries.push_back({row, column, mv});
                if (check_dense_reference) {
                    const size_t index = static_cast<size_t>(row) * static_cast<size_t>(n) +
                                         static_cast<size_t>(column);
                    dense_reference_k[index] += kv;
                    dense_reference_m[index] += mv;
                }
            }
        }
    }

    const size_t triplet_capacity_bytes =
        k_entries.capacity() * sizeof(Triplet) + m_entries.capacity() * sizeof(Triplet);
    const SparseMatrix Ks = SparseMatrix::from_triplets(n, std::move(k_entries));
    const SparseMatrix Ms = SparseMatrix::from_triplets(n, std::move(m_entries));
    result.matrix_audit.full_dof = static_cast<size_t>(ndof_full);
    result.matrix_audit.reduced_dof = static_cast<size_t>(n);
    result.matrix_audit.assembly_dense_bytes = 0u;
    result.matrix_audit.assembly_sparse_bytes = Ks.storage_bytes() + Ms.storage_bytes();
    result.matrix_audit.assembly_triplet_capacity_bytes = triplet_capacity_bytes;
    result.matrix_audit.stiffness = audit_matrix(Ks);
    result.matrix_audit.mass = audit_matrix(Ms);
    // Diagnostic common RCM ordering: one permutation must be applicable to
    // both K and M in the generalized eigenproblem. No production solve uses
    // this ordering yet.
    const auto common_graph = detail::structural_union(Ks, Ms);
    const auto rcm_order = detail::reverse_cuthill_mckee_order(common_graph);
    const auto K_rcm = detail::symmetric_permute(Ks, rcm_order);
    const auto M_rcm = detail::symmetric_permute(Ms, rcm_order);
    const auto K_rcm_metrics = detail::structural_metrics(K_rcm);
    const auto M_rcm_metrics = detail::structural_metrics(M_rcm);
    result.matrix_audit.stiffness.rcm_max_bandwidth = K_rcm_metrics.max_bandwidth;
    result.matrix_audit.stiffness.rcm_profile = K_rcm_metrics.profile;
    result.matrix_audit.mass.rcm_max_bandwidth = M_rcm_metrics.max_bandwidth;
    result.matrix_audit.mass.rcm_profile = M_rcm_metrics.profile;
    result.matrix_audit.available = true;
    result.matrix_audit.dense_reference_checked = check_dense_reference;
    if (check_dense_reference) {
        double k_error = 0.0, m_error = 0.0;
        double k_scale = 0.0, m_scale = 0.0;
        for (int row = 0; row < n; ++row) {
            for (int column = 0; column < n; ++column) {
                const size_t index = static_cast<size_t>(row) * static_cast<size_t>(n) +
                                     static_cast<size_t>(column);
                const double kv = Ks.get(row, column), mv = Ms.get(row, column);
                k_error = std::max(k_error, std::abs(kv - dense_reference_k[index]));
                m_error = std::max(m_error, std::abs(mv - dense_reference_m[index]));
                k_scale = std::max(k_scale, std::abs(dense_reference_k[index]));
                m_scale = std::max(m_scale, std::abs(dense_reference_m[index]));
            }
        }
        result.matrix_audit.stiffness_dense_reference_max_relative_error =
            k_error / std::max(k_scale, std::numeric_limits<double>::min());
        result.matrix_audit.mass_dense_reference_max_relative_error =
            m_error / std::max(m_scale, std::numeric_limits<double>::min());
    }
    if (!quiet_) {
        const auto print_stats = [](const char* name, const MatrixAuditStats& stats) {
            std::cerr << "Matrix audit " << name << ": n=" << stats.dimension
                      << " nnz=" << stats.exact_nonzeros
                      << " density=" << stats.density
                      << " dense_bytes=" << stats.dense_bytes
                      << " csr_bytes=" << stats.estimated_csr_bytes
                      << " relative_asymmetry=" << stats.relative_asymmetry_frobenius
                      << " diag_min/max/range=" << stats.positive_diagonal_min << '/'
                      << stats.positive_diagonal_max << '/' << stats.positive_diagonal_dynamic_range
                      << " nonpositive_diag=" << stats.nonpositive_diagonal_count
                      << " candidate_pairs<=1e-12/1e-10/1e-8="
                      << stats.candidate_pairs_below_1e12 << '/'
                      << stats.candidate_pairs_below_1e10 << '/'
                      << stats.candidate_pairs_below_1e8
                      << " structural_pairs=" << stats.symmetric_offdiagonal_pairs
                      << " structural_asymmetry=" << stats.structural_asymmetry_pairs
                      << " max_bandwidth=" << stats.max_bandwidth
                      << " profile=" << stats.profile
                      << " mean_row_nnz=" << stats.mean_row_nnz << '\n';
        };
        std::cerr << "Matrix audit full_dof=" << result.matrix_audit.full_dof
                  << " reduced_dof=" << result.matrix_audit.reduced_dof
                  << " assembly_dense_bytes=0 assembly_sparse_bytes="
                  << result.matrix_audit.assembly_sparse_bytes
                  << " triplet_capacity_bytes=" << result.matrix_audit.assembly_triplet_capacity_bytes << '\n';
        print_stats("K", result.matrix_audit.stiffness);
        print_stats("M", result.matrix_audit.mass);
    }

    // Dense Jacobi remains the reference for small systems. Larger systems use
    // dependency-free sparse shifted-inverse subspace iteration; no n-by-n dense
    // matrices are materialized on that path.
    const bool use_sparse_iteration = n > sparse_solver_threshold_;
    Matrix L;
    std::vector<double> eval;
    Matrix evec;
    std::vector<std::vector<double>> sparse_mode_vectors;
    if (use_sparse_iteration) {
        std::vector<Triplet> k_scaled_entries, m_scaled_entries;
        k_scaled_entries.reserve(Ks.values.size());
        m_scaled_entries.reserve(Ms.values.size());
        for (int row = 0; row < n; ++row) {
            for (size_t p = Ks.row_offsets[static_cast<size_t>(row)];
                 p < Ks.row_offsets[static_cast<size_t>(row) + 1u]; ++p)
                k_scaled_entries.push_back({row, Ks.columns[p], Ks.values[p] / material.E});
            for (size_t p = Ms.row_offsets[static_cast<size_t>(row)];
                 p < Ms.row_offsets[static_cast<size_t>(row) + 1u]; ++p)
                m_scaled_entries.push_back({row, Ms.columns[p], Ms.values[p] / material.rho});
        }
        const SparseMatrix Kscaled = SparseMatrix::from_triplets(n, std::move(k_scaled_entries));
        const SparseMatrix Mscaled = SparseMatrix::from_triplets(n, std::move(m_scaled_entries));
        const int wanted = std::min(n, requested_modes_ + 8);
        auto iterative = detail::solve_lowest_sparse_modes(Kscaled, Mscaled, wanted, 500, 1e-5, use_rcm_ordering_);
        bool used_rcm = use_rcm_ordering_;
        // RCM changes the numerical IC(0) factor even though it preserves the
        // physical matrix. It can substantially accelerate PCG, but it is not
        // guaranteed to converge for every spectrum. If the RCM path fails its
        // own convergence/backward-error contract, retry in the original order
        // rather than weakening the accuracy requirement or changing physics.
        if (!iterative.converged && use_rcm_ordering_) {
            iterative = detail::solve_lowest_sparse_modes(Kscaled, Mscaled, wanted, 500, 1e-5, false);
            used_rcm = false;
        }
        if (!iterative.converged) {
            if (!quiet_) std::cerr << "Internal FEM: sparse eigensolver did not converge (iterations="
                               << iterative.iterations << ", returned=" << iterative.eigenvalues.size()
                               << ", relative_residual=" << iterative.final_relative_residual
                               << ", eigenvalue_change=" << iterative.final_eigenvalue_change << ")\n";
            return false;
        }
        eval = iterative.eigenvalues;
        sparse_mode_vectors = iterative.vectors;
        result.matrix_audit.production_rcm_ordering = used_rcm;
        result.matrix_audit.stiffness.ic0_factor_nonzeros = iterative.ic0_nonzeros;
        result.matrix_audit.sparse_solver_iterations = iterative.iterations;
        result.matrix_audit.sparse_pcg_iterations = iterative.pcg_iterations;
        if (!quiet_) std::cerr << "Internal FEM sparse eigensolver: iterations="
                               << iterative.iterations << " returned=" << eval.size() << '\n';
    } else {
        Matrix Kr = sparse_to_dense(Ks);
        Matrix Mr = sparse_to_dense(Ms);
        Matrix Kr_scaled = Kr;
        Matrix Mr_scaled = Mr;
        for (double& x : Kr_scaled.a) x /= material.E;
        for (double& x : Mr_scaled.a) x /= material.rho;
        if (!cholesky(Mr_scaled,L)) {
            if (!quiet_) std::cerr << "Internal FEM: mass matrix is not positive definite\n";
            return false;
        }
        Matrix A = generalized_to_standard(Kr_scaled,L);
        if (!jacobi(A,eval,evec)) return false;
    }

        std::vector<double> q_radial(static_cast<size_t>(n), 0.0);
    std::vector<double> q_axial(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        const int full_dof = free_dof[static_cast<size_t>(i)];
        const bool radial = (full_dof % 2) == 0;
        if (radial) q_radial[static_cast<size_t>(i)] = 1.0;
        else q_axial[static_cast<size_t>(i)] = 1.0;
    }

    // Excitation is supplied through a neutral influence contract. FEM owns
    // the mass projection; it does not construct a strike/rotation field.
    if (problem.influence_provider == nullptr) {
        if (!quiet_) std::cerr << "Internal FEM: a modal influence provider is required\n";
        return false;
    }
    const ModalInfluence influence = problem.influence_provider->make(mesh);
    if (influence.radial.size() != mesh.nodes.size() ||
        influence.axial.size() != mesh.nodes.size()) {
        if (!quiet_) std::cerr << "Internal FEM: modal influence has invalid nodal dimensions\n";
        return false;
    }
    std::vector<double> q_influence(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        const int full_dof = free_dof[static_cast<size_t>(i)];
        const int node_index = full_dof / 2;
        q_influence[static_cast<size_t>(i)] =
            (full_dof % 2 == 0)
                ? influence.radial[static_cast<size_t>(node_index)]
                : influence.axial[static_cast<size_t>(node_index)];
    }



    double spectral_scale = 0.0;
    for (double value : eval) spectral_scale = std::max(spectral_scale, std::abs(value));
    const double zero_eigenvalue_tolerance = std::max(
        spectral_scale * 1e-12, std::numeric_limits<double>::min());

    std::vector<std::vector<double>> accepted_modes;
    for (int i = 0; i < static_cast<int>(eval.size()) && static_cast<int>(result.modes.size()) < requested_modes_; ++i) {
        // Classify numerical rigid/zero modes relative to the spectrum, not by an
        // absolute frequency. A fixed 0.5-Hz cutoff incorrectly discards valid
        // low-frequency modes when material stiffness or geometry is rescaled.
        if (!(eval[i] > zero_eigenvalue_tolerance) || !std::isfinite(eval[i])) continue;
        const double lambda = eval[i] * material.E / material.rho;
        if (!(lambda > 0.0) || !std::isfinite(lambda)) continue;
        const double freq = std::sqrt(lambda) / (2.0 * PI);

        std::vector<double> phi;
        if (use_sparse_iteration) {
            if (static_cast<size_t>(i) >= sparse_mode_vectors.size()) continue;
            phi = sparse_mode_vectors[static_cast<size_t>(i)];
        } else {
            std::vector<double> y(static_cast<size_t>(n));
            for (int k = 0; k < n; ++k) y[static_cast<size_t>(k)] = evec(k, i);
            solve_upper_transpose(L, y, phi);
        }

        // The eigensolver uses M/rho. Restore physical mass normalization
        // so participation factors and effective masses have a well-defined
        // convention independent of the dimensionless reduction.
        const auto physical_mphi = Ms.multiply(phi);
        double physical_modal_mass = 0.0;
        for (int row = 0; row < n; ++row)
            physical_modal_mass += phi[static_cast<size_t>(row)] * physical_mphi[static_cast<size_t>(row)];
        if (!(physical_modal_mass > 0.0) || !std::isfinite(physical_modal_mass)) continue;
        const double norm = 1.0 / std::sqrt(physical_modal_mass);
        for (double& value : phi) value *= norm;
        const auto normalized_mphi = Ms.multiply(phi);
        double normalized_modal_mass = 0.0;
        for (int row = 0; row < n; ++row)
            normalized_modal_mass += phi[static_cast<size_t>(row)] * normalized_mphi[static_cast<size_t>(row)];

        auto participation = [&](const std::vector<double>& q) {
            double value = 0.0;
            const auto mphi = Ms.multiply(phi);
            for (int row = 0; row < n; ++row) value += q[static_cast<size_t>(row)] * mphi[static_cast<size_t>(row)];
            return value;
        };

        const auto& profile = problem.mesh.outer_profile;
        const size_t outer_count = profile.size();
        std::vector<int> full_to_free(static_cast<size_t>(ndof_full), -1);
        for (int ri = 0; ri < n; ++ri) full_to_free[static_cast<size_t>(free_dof[ri])] = ri;

        // Validate the physical generalized eigenpair after normalization.
        // This is a diagnostic, not a substitute for independent mesh convergence.
        const auto kphi = Ks.multiply(phi);
        const auto mphi = Ms.multiply(phi);
        double k_norm_sq = 0.0, m_norm_sq = 0.0, residual_norm_sq = 0.0;
        for (int row = 0; row < n; ++row) {
            const double residual = kphi[static_cast<size_t>(row)] -
                                    lambda * mphi[static_cast<size_t>(row)];
            k_norm_sq += kphi[static_cast<size_t>(row)] * kphi[static_cast<size_t>(row)];
            m_norm_sq += mphi[static_cast<size_t>(row)] * mphi[static_cast<size_t>(row)];
            residual_norm_sq += residual * residual;
        }
        const double residual_denominator = std::sqrt(k_norm_sq) +
                                            std::abs(lambda) * std::sqrt(m_norm_sq);
        const double eigen_residual_relative = std::sqrt(residual_norm_sq) /
            std::max(residual_denominator, std::numeric_limits<double>::min());
        double mass_orthogonality_error = 0.0;
        for (const auto& previous : accepted_modes) {
            const auto m_previous = Ms.multiply(previous);
            double overlap = 0.0;
            for (int row = 0; row < n; ++row)
                overlap += phi[static_cast<size_t>(row)] * m_previous[static_cast<size_t>(row)];
            mass_orthogonality_error = std::max(mass_orthogonality_error, std::abs(overlap));
        }

        ModalMode mode;
        mode.eigen_residual_relative = eigen_residual_relative;
        mode.modal_mass = normalized_modal_mass;
        mode.mass_orthogonality_error = mass_orthogonality_error;
        accepted_modes.push_back(phi);
        mode.mode_number = static_cast<int>(result.modes.size()) + 2;
        mode.eigenvalue = lambda;
        mode.frequency_hz = freq;
        mode.part_x = participation(q_radial);
        mode.part_y = participation(q_axial);
        mode.part_z = 0.0;
        mode.part_rx = 0.0;
        mode.part_ry = 0.0;
        mode.part_rz = participation(q_influence);
        mode.mass_x = mode.part_x * mode.part_x;
        mode.mass_y = mode.part_y * mode.part_y;
        mode.mass_z = 0.0;
        mode.mass_rx = 0.0;
        mode.mass_ry = 0.0;
        // This squared generalized influence is a diagnostic normalization, not a
        // physical 3-D rotational inertia. The meaning is defined by the
        // selected influence provider.
        mode.mass_rz = mode.part_rz * mode.part_rz;
        mode.effective_mass_norm = std::abs(mode.mass_rz);

        // Export only geometry-independent boundary samples when explicitly
        // requested. The FEM solver does not construct strike normals or any
        // other excitation-specific quantity. Excitation modules may project
        // arbitrary loads onto these modal boundary displacements.
        if (capture_mode_shapes_) {
            mode.outer_surface_radial.assign(profile.size(), 0.0);
            mode.outer_surface_axial.assign(profile.size(), 0.0);
            for (size_t k = 0; k < profile.size(); ++k) {
                const int rr = full_to_free[2 * k];
                const int zz = full_to_free[2 * k + 1];
                mode.outer_surface_radial[k] = (rr >= 0) ? phi[static_cast<size_t>(rr)] : 0.0;
                mode.outer_surface_axial[k] = (zz >= 0) ? phi[static_cast<size_t>(zz)] : 0.0;
            }
        }

        // Acoustic radiation coefficient for the axisymmetric outer wall.
        // The eigenvector is mass-normalized, so this is the displaced volume
        // per unit modal coordinate: Qd = integral_S (phi . n_out) dS.
        // For the profile orientation z increasing from the axis to the rim,
        // n_out dS = (dz, -dr). This lets us integrate the surface flux
        // directly, without introducing a secant-normal approximation.
        // The FE displacement is linear on each profile segment, so the
        // following segment integral is exact for that piecewise-linear field.
        if (outer_count >= 2 && outer_count <= mesh.nodes.size()) {
            double qd = 0.0;
            double dz_moment = 0.0;
            bool valid_qd = true;
            // Use the mouth plane as the acoustic origin. This is the natural
            // reference for the open bell and also keeps the first moment
            // numerically well-scaled. For Qd ~= 0 the result is origin
            // independent to numerical precision.
            const double z_ref = profile.back().z;
            for (size_t k = 0; k + 1 < outer_count; ++k) {
                const Point2D& p0 = profile[k];
                const Point2D& p1 = profile[k + 1];
                const double dr = p1.r - p0.r;
                const double dz = p1.z - p0.z;
                if (!std::isfinite(dr) || !std::isfinite(dz) ||
                    !(std::hypot(dr, dz) > 0.0)) {
                    valid_qd = false;
                    break;
                }

                const size_t node0 = k;
                const size_t node1 = k + 1;
                const int r0 = full_to_free[2 * node0];
                const int z0 = full_to_free[2 * node0 + 1];
                const int r1 = full_to_free[2 * node1];
                const int z1 = full_to_free[2 * node1 + 1];
                if (r0 < 0 || z0 < 0 || r1 < 0 || z1 < 0) {
                    // The radial DOF at the symmetry axis is constrained; its
                    // contribution is exactly zero by axis regularity.
                    if (k != 0 || z0 < 0 || r1 < 0 || z1 < 0) {
                        valid_qd = false;
                        break;
                    }
                }

                const double ur0 = (r0 >= 0) ? phi[static_cast<size_t>(r0)] : 0.0;
                const double uz0 = (z0 >= 0) ? phi[static_cast<size_t>(z0)] : 0.0;
                const double ur1 = (r1 >= 0) ? phi[static_cast<size_t>(r1)] : 0.0;
                const double uz1 = (z1 >= 0) ? phi[static_cast<size_t>(z1)] : 0.0;

                // With t in [0,1], r=r0+t*dr and u=u0+t*(u1-u0).
                // Let F(t)=r(t)*(ur(t)*dz-uz(t)*dr)=a0+a1*t+a2*t^2.
                // Qd is 2*pi*integral F dt.
                const double a0 = p0.r * (ur0 * dz - uz0 * dr);
                const double a1 = dr * (ur0 * dz - uz0 * dr)
                                + p0.r * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
                const double a2 = dr * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
                const double i0 = a0 + 0.5 * a1 + a2 / 3.0;
                qd += 2.0 * PI * i0;

                // First axial moment: Dz = 2*pi*integral z(t)*F(t) dt.
                // z(t)=z0+t*dz, so the integral is exact cubic integration.
                const double i1 = p0.z * i0
                                + dz * (0.5 * a0 + a1 / 3.0 + a2 / 4.0)
                                - z_ref * i0;
                dz_moment += 2.0 * PI * i1;
            }
            if (valid_qd && std::isfinite(qd) && std::isfinite(dz_moment)) {
                mode.acoustic_volume_displacement = qd;
                mode.acoustic_dipole_moment_z = dz_moment;
            }
        }

        if (capture_mode_shapes_) {
            mode.nodal_radial.assign(mesh.nodes.size(), 0.0);
            mode.nodal_axial.assign(mesh.nodes.size(), 0.0);
            for (size_t node_index = 0; node_index < mesh.nodes.size(); ++node_index) {
                const int radial_dof = static_cast<int>(2 * node_index);
                const int axial_dof = radial_dof + 1;
                const int radial_free = full_to_free[static_cast<size_t>(radial_dof)];
                const int axial_free = full_to_free[static_cast<size_t>(axial_dof)];
                if (radial_free >= 0) mode.nodal_radial[node_index] = phi[static_cast<size_t>(radial_free)];
                if (axial_free >= 0) mode.nodal_axial[node_index] = phi[static_cast<size_t>(axial_free)];
            }
        }

        result.modes.push_back(mode);
    }

    double max_mechanical_excitation = 0.0;
    double max_acoustic_coefficient = 0.0;
    for (const auto& mode : result.modes) {
        max_mechanical_excitation = std::max(
            max_mechanical_excitation,
            std::sqrt(std::max(0.0, mode.effective_mass_norm)));
        const double omega = 2.0 * PI * mode.frequency_hz;
        const double coefficient =
            -mode.part_rz * mode.acoustic_volume_displacement * omega;
        if (std::isfinite(coefficient))
            max_acoustic_coefficient = std::max(max_acoustic_coefficient,
                                                std::abs(coefficient));
    }
    for (auto& mode : result.modes) {
        const double excitation = std::sqrt(std::max(0.0, mode.effective_mass_norm));
        mode.modal_excitation = excitation;

        const double omega = 2.0 * PI * mode.frequency_hz;
        const double coefficient =
            -mode.part_rz * mode.acoustic_volume_displacement * omega;
        if (std::isfinite(coefficient)) {
            mode.acoustic_pressure_coefficient = coefficient;
            mode.acoustic_amplitude = (max_acoustic_coefficient > 1e-14)
                ? std::abs(coefficient) / max_acoustic_coefficient
                : 0.0;
        }

        // Production amplitude is acoustic when the radiation integral is
        // available. Mechanical participation remains available separately
        // in modal_excitation/effective_mass_norm.
        if (max_acoustic_coefficient > 1e-14) {
            mode.amplitude = mode.acoustic_amplitude;
        } else {
            mode.amplitude = (max_mechanical_excitation > 1e-10)
                ? excitation / max_mechanical_excitation
                : 1.0;
        }
    }

    result.sort_by_frequency();
    if(!result.modes.empty()) result.fundamental_frequency=result.modes.front().frequency_hz;
    if(!quiet_) std::cout<<"Internal FEM modes: "<<result.modes.size()<<"\n";
    return !result.modes.empty();
}

} // namespace BellOptimizer
