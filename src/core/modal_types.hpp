#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace BellOptimizer {

struct MatrixAuditStats {
    size_t dimension = 0;
    size_t dense_bytes = 0;
    size_t exact_nonzeros = 0;
    size_t estimated_csr_bytes = 0;
    double density = 0.0;
    double relative_asymmetry_frobenius = 0.0;
    // Diagonal scaling diagnostic. A large ratio can indicate poor scaling or
    // near-mechanisms, but is not by itself a proof of ill-conditioning.
    double positive_diagonal_min = 0.0;
    double positive_diagonal_max = 0.0;
    double positive_diagonal_dynamic_range = 0.0;
    size_t nonpositive_diagonal_count = 0;
    // Symmetry-preserving candidate pair counts; diagnostics only, no entries removed.
    size_t candidate_pairs_below_1e12 = 0;
    size_t candidate_pairs_below_1e10 = 0;
    size_t candidate_pairs_below_1e8 = 0;
    // Structural graph diagnostics. These describe the FEM connectivity, not
    // numerical coefficient magnitudes, and are useful for sparse ordering.
    size_t symmetric_offdiagonal_pairs = 0;
    size_t structural_asymmetry_pairs = 0;
    size_t max_bandwidth = 0;
    size_t profile = 0;
    double mean_row_nnz = 0.0;
    // Diagnostic ordering metrics; production sparse FEM uses the configured RCM ordering.
    size_t rcm_max_bandwidth = 0;
    size_t rcm_profile = 0;
    // IC(0) factor nonzeros for the production sparse solve. Useful for
    // ordering diagnostics; zero when the sparse solver was not used.
    size_t ic0_factor_nonzeros = 0;
};

struct FEMMatrixAudit {
    size_t full_dof = 0;
    size_t reduced_dof = 0;
    size_t assembly_dense_bytes = 0; // zero when global assembly is sparse
    size_t assembly_sparse_bytes = 0;
    size_t assembly_triplet_capacity_bytes = 0; // transient COO workspace, before CSR compression
    bool dense_reference_checked = false;
    double stiffness_dense_reference_max_relative_error = 0.0;
    double mass_dense_reference_max_relative_error = 0.0;
    MatrixAuditStats stiffness;
    MatrixAuditStats mass;
    bool available = false;
    bool production_rcm_ordering = false;
    int sparse_solver_iterations = 0;
    std::size_t sparse_pcg_iterations = 0;
};

// Neutral generalized modal influence field. It describes a displacement/load
// influence in the reduced FEM DOF ordering; it does not prescribe a physical
// excitation mechanism. Concrete strike/rotation models live outside FEM.
struct ModalInfluence {
    std::vector<double> radial;
    std::vector<double> axial;
    std::string id;
    double reference_z = 0.0;
};

struct ModalMode {
    int mode_number = 0;
    double eigenvalue = 0.0;
    double frequency_hz = 0.0;
    // Normalized audible modal amplitude. When the acoustic surface coefficient
    // is available this is |C_m| normalized across returned modes, where
    // C_m = -part_rz * Qd * omega for the current excitation convention.
    // Legacy mechanical fallback is used only when acoustic data is unavailable.
    double amplitude = 0.0;
    double decay_time = 0.0;
    // Numerical validation diagnostics from the internal FEM eigensolver.
    // Residual is ||K phi - lambda M phi|| / (||K phi|| + |lambda| ||M phi||).
    double eigen_residual_relative = 0.0;
    // Physical modal mass phi^T M phi; should be 1 for internal mass-normalized modes.
    double modal_mass = 0.0;
    // Maximum |phi_i^T M phi_j| against previously accepted modes in this solve.
    double mass_orthogonality_error = 0.0;

    double part_x = 0.0, part_y = 0.0, part_z = 0.0;
    double part_rx = 0.0, part_ry = 0.0, part_rz = 0.0;

    double mass_x = 0.0, mass_y = 0.0, mass_z = 0.0;
    double mass_rx = 0.0, mass_ry = 0.0, mass_rz = 0.0;
    // Legacy squared RZ/meridional-proxy participation, retained for file/API
    // compatibility. Do not assume it is a physical rotational effective
    // inertia for the axisymmetric model; its value depends on the reference
    // origin and the excitation field is not a 3-D axisymmetric rigid rotation.
    double effective_mass_norm = 0.0;
    // Surface volume-displacement coefficient of the mass-normalized mode:
    // Qd = integral_S (phi . n) dS. For the axisymmetric outer wall
    // this is the volume displaced per unit modal coordinate.
    double acoustic_volume_displacement = 0.0;
    // First axial acoustic radiation moment about the selected acoustic origin:
    // Dz = integral_S (z - z_ref) (phi . n) dS.
    double acoustic_dipole_moment_z = 0.0;
    // Signed far-field monopole pressure coefficient per unit common excitation
    // impulse. Common factors (air density, distance, excitation magnitude)
    // cancel in modal ratios but the sign is retained for phase.
    double acoustic_pressure_coefficient = 0.0;
    // Absolute value of the acoustic pressure coefficient, normalized across
    // returned modes. This is the quantity used by the production objective.
    double acoustic_amplitude = 0.0;
    // Mechanical modal excitation proxy derived from the squared generalized
    // influence selected for the solve. It is separate from acoustics.
    double modal_excitation = 0.0;
    // Mass-normalized modal displacement samples on the outer meridional
    // boundary. These are boundary data, not an excitation model. They allow
    // independent excitation modules to project arbitrary loads onto the
    // FEM mode shapes without embedding strike assumptions in the solver.
    // Values may change sign with the arbitrary eigenvector sign.
    std::vector<double> outer_surface_radial;
    std::vector<double> outer_surface_axial;
    // Optional complete Fourier boundary representation from direct vector FEM.
    // Layout is [profile_point][harmonic_basis]; each component array has equal
    // length profile_point_count * outer_surface_harmonics.size(). Positive
    // labels denote cosine terms; negative labels denote sine terms. The legacy
    // radial/axial arrays above remain the axisymmetric (basis-0) view.
    std::vector<int> outer_surface_harmonics;
    std::vector<double> outer_surface_fourier_radial;
    std::vector<double> outer_surface_fourier_azimuthal;
    std::vector<double> outer_surface_fourier_axial;
    // Full-mesh eigenvector samples in mesh.nodes order, used for independent FEM validation.
    std::vector<double> nodal_radial;
    std::vector<double> nodal_axial;
    // Complete direct-vector-FEM eigenfield in global DOF order:
    // [node][Fourier basis][ur, uphi, uz], mass-normalized. The associated
    // diagonal row-sum weights are a lumped-M approximation for diagnostics;
    // they are not a replacement for the consistent FEM mass matrix.
    std::vector<double> fem_displacement;
    std::vector<double> fem_mass_lumped_weights;
    std::vector<double> fem_node_r;
    std::vector<double> fem_node_z;
    std::vector<int> fem_harmonics;
};

struct ModalData {
    std::vector<ModalMode> modes;
    FEMMatrixAudit matrix_audit;
    double fundamental_frequency = 0.0;
    std::string source_file;

    double sample_rate = 48000.0;
    double duration = 3.0;
    double t60_base = 2.0;
    double target_rms = 0.70;
    std::string material;

    void clear() { modes.clear(); matrix_audit = {}; fundamental_frequency = 0.0; source_file.clear(); }

    void sort_by_frequency() {
        std::sort(modes.begin(), modes.end(),
                  [](const ModalMode& a, const ModalMode& b) {
                      return a.frequency_hz < b.frequency_hz;
                  });
    }

    void normalize_amplitudes(double max_amp = 1.0) {
        if (modes.empty()) return;
        double max_amp_found = 0.0;
        for (const auto& mode : modes) max_amp_found = std::max(max_amp_found, mode.amplitude);
        if (max_amp_found > 1e-10) {
            const double scale = max_amp / max_amp_found;
            for (auto& mode : modes) mode.amplitude *= scale;
        }
    }

    size_t size() const { return modes.size(); }
    bool empty() const { return modes.empty(); }
};


} // namespace BellOptimizer
