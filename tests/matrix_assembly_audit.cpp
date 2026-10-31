#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <cmath>
#include <iostream>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 2;

    BellOptimizer::InternalModalSolver solver(false, 6);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "matrix_assembly_audit", &modal_influence_1};
    BellOptimizer::ModalData result;
    if (!solver.solve(problem, result)) return 3;
    const auto& audit = result.matrix_audit;
    if (!audit.available || audit.full_dof == 0 || audit.reduced_dof == 0) return 4;
    if (audit.stiffness.dimension != audit.reduced_dof || audit.mass.dimension != audit.reduced_dof) return 5;
    if (audit.stiffness.exact_nonzeros == 0 || audit.mass.exact_nonzeros == 0) return 6;
    if (!(audit.stiffness.density > 0.0 && audit.stiffness.density <= 1.0) ||
        !(audit.mass.density > 0.0 && audit.mass.density <= 1.0)) return 7;
    if (audit.stiffness.relative_asymmetry_frobenius > 1e-12 ||
        audit.mass.relative_asymmetry_frobenius > 1e-12) return 8;
    if (audit.stiffness.estimated_csr_bytes >= audit.stiffness.dense_bytes ||
        audit.mass.estimated_csr_bytes >= audit.mass.dense_bytes) return 9;
    if (audit.assembly_dense_bytes != 0u || audit.assembly_sparse_bytes == 0u ||
        audit.assembly_triplet_capacity_bytes == 0u) return 10;
    if (audit.assembly_sparse_bytes >= 2u * audit.reduced_dof * audit.reduced_dof * sizeof(double)) return 11;
    if (!audit.dense_reference_checked ||
        audit.stiffness_dense_reference_max_relative_error > 1e-12 ||
        audit.mass_dense_reference_max_relative_error > 1e-12) return 12;
    if (audit.stiffness.structural_asymmetry_pairs != 0 ||
        audit.mass.structural_asymmetry_pairs != 0) return 13;
    if (audit.stiffness.symmetric_offdiagonal_pairs == 0 ||
        audit.mass.symmetric_offdiagonal_pairs == 0) return 14;
    if (audit.stiffness.max_bandwidth == 0 || audit.mass.max_bandwidth == 0) return 15;

    std::cout << "matrix_assembly_audit=PASS\n"
              << "full_dof=" << audit.full_dof << '\n'
              << "reduced_dof=" << audit.reduced_dof << '\n'
              << "assembly_dense_bytes=" << audit.assembly_dense_bytes << '\n'
              << "assembly_sparse_bytes=" << audit.assembly_sparse_bytes << '\n'
              << "assembly_triplet_capacity_bytes=" << audit.assembly_triplet_capacity_bytes << '\n'
              << "K_density=" << audit.stiffness.density << '\n'
              << "M_density=" << audit.mass.density << '\n'
              << "K_dense_bytes=" << audit.stiffness.dense_bytes << '\n'
              << "K_csr_estimated_bytes=" << audit.stiffness.estimated_csr_bytes << '\n'
              << "M_dense_bytes=" << audit.mass.dense_bytes << '\n'
              << "M_csr_estimated_bytes=" << audit.mass.estimated_csr_bytes << '\n'
              << "K_relative_asymmetry=" << audit.stiffness.relative_asymmetry_frobenius << '\n'
              << "M_relative_asymmetry=" << audit.mass.relative_asymmetry_frobenius << '\n'
              << "dense_reference_checked=" << audit.dense_reference_checked << '\n'
              << "K_dense_reference_max_relative_error=" << audit.stiffness_dense_reference_max_relative_error << '\n'
              << "M_dense_reference_max_relative_error=" << audit.mass_dense_reference_max_relative_error << '\n'
              << "K_candidate_pairs_below_1e-12=" << audit.stiffness.candidate_pairs_below_1e12 << '\n'
              << "K_candidate_pairs_below_1e-10=" << audit.stiffness.candidate_pairs_below_1e10 << '\n'
              << "K_candidate_pairs_below_1e-8=" << audit.stiffness.candidate_pairs_below_1e8 << '\n'
              << "M_candidate_pairs_below_1e-12=" << audit.mass.candidate_pairs_below_1e12 << '\n'
              << "M_candidate_pairs_below_1e-10=" << audit.mass.candidate_pairs_below_1e10 << '\n'
              << "M_candidate_pairs_below_1e-8=" << audit.mass.candidate_pairs_below_1e8 << '\n'
              << "K_symmetric_offdiagonal_pairs=" << audit.stiffness.symmetric_offdiagonal_pairs << '\n'
              << "K_structural_asymmetry_pairs=" << audit.stiffness.structural_asymmetry_pairs << '\n'
              << "K_max_bandwidth=" << audit.stiffness.max_bandwidth << '\n'
              << "K_profile=" << audit.stiffness.profile << '\n'
              << "K_mean_row_nnz=" << audit.stiffness.mean_row_nnz << '\n'
              << "M_symmetric_offdiagonal_pairs=" << audit.mass.symmetric_offdiagonal_pairs << '\n'
              << "M_structural_asymmetry_pairs=" << audit.mass.structural_asymmetry_pairs << '\n'
              << "M_max_bandwidth=" << audit.mass.max_bandwidth << '\n'
              << "M_profile=" << audit.mass.profile << '\n'
              << "M_mean_row_nnz=" << audit.mass.mean_row_nnz << '\n'
              << "M_rcm_max_bandwidth=" << audit.mass.rcm_max_bandwidth << '\n'
              << "M_rcm_profile=" << audit.mass.rcm_profile << '\n';
    return 0;
}
