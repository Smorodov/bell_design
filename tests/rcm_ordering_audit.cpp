#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "fem/sparse_matrix.hpp"

#include <cmath>
#include <iostream>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;

static bool same_dense_under_permutation(const SparseMatrix& a, const SparseMatrix& b,
                                         const std::vector<int>& order, double tol) {
    if (a.n != b.n || static_cast<int>(order.size()) != a.n) return false;
    for (int i = 0; i < a.n; ++i)
        for (int j = 0; j < a.n; ++j) {
            const double expected = a.get(order[static_cast<size_t>(i)], order[static_cast<size_t>(j)]);
            if (std::abs(expected - b.get(i, j)) > tol * std::max(1.0, std::abs(expected))) return false;
        }
    return true;
}

int main() {
    // A deliberately ordered path with a long-distance edge. RCM must reduce
    // the structural bandwidth while preserving every coefficient.
    std::vector<SparseTriplet> entries;
    const int n = 12;
    for (int i = 0; i < n; ++i) entries.push_back({i, i, 2.0});
    for (int i = 0; i + 1 < n; ++i) {
        entries.push_back({i, i + 1, -1.0});
        entries.push_back({i + 1, i, -1.0});
    }
    entries.push_back({0, 6, -0.25}); entries.push_back({6, 0, -0.25});
    entries.push_back({6, 11, -0.25}); entries.push_back({11, 6, -0.25});
    const SparseMatrix A = SparseMatrix::from_triplets(n, entries);
    const auto before = BellOptimizer::detail::structural_metrics(A);
    const auto order = BellOptimizer::detail::reverse_cuthill_mckee_order(A);
    const auto reordered = BellOptimizer::detail::symmetric_permute(A, order);
    const auto after = BellOptimizer::detail::structural_metrics(reordered);
    if (static_cast<int>(order.size()) != n) return 2;
    if (!same_dense_under_permutation(A, reordered, order, 1e-14)) return 3;
    if (reordered.values.size() != A.values.size()) return 4;
    if (after.max_bandwidth > before.max_bandwidth) return 5;

    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 31;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 6;
    BellOptimizer::InternalModalSolver solver(true, 6);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "rcm_ordering_audit", &modal_influence_1};
    BellOptimizer::ModalData result;
    if (!solver.solve(problem, result)) return 7;
    const auto& audit = result.matrix_audit;
    if (!audit.available) return 8;
    if (audit.stiffness.rcm_max_bandwidth > audit.stiffness.max_bandwidth ||
        audit.mass.rcm_max_bandwidth > audit.mass.max_bandwidth) return 9;
    if (audit.stiffness.rcm_profile > audit.stiffness.profile ||
        audit.mass.rcm_profile > audit.mass.profile) return 10;

    std::cout << "rcm_ordering_audit=PASS\n"
              << "synthetic_bandwidth_before=" << before.max_bandwidth << '\n'
              << "synthetic_bandwidth_after=" << after.max_bandwidth << '\n'
              << "synthetic_profile_before=" << before.profile << '\n'
              << "synthetic_profile_after=" << after.profile << '\n'
              << "K_bandwidth_before=" << audit.stiffness.max_bandwidth << '\n'
              << "K_bandwidth_after=" << audit.stiffness.rcm_max_bandwidth << '\n'
              << "K_profile_before=" << audit.stiffness.profile << '\n'
              << "K_profile_after=" << audit.stiffness.rcm_profile << '\n'
              << "M_bandwidth_before=" << audit.mass.max_bandwidth << '\n'
              << "M_bandwidth_after=" << audit.mass.rcm_max_bandwidth << '\n'
              << "M_profile_before=" << audit.mass.profile << '\n'
              << "M_profile_after=" << audit.mass.rcm_profile << '\n';
    return 0;
}
