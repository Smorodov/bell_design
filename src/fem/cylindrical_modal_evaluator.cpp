#include "cylindrical_modal_evaluator.hpp"
#include "cylindrical_global_modes.hpp"
#include "cylindrical_modal_tracking.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BellOptimizer::FEM {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;

bool finite_positive(double v) { return std::isfinite(v) && v > 0.0; }

} // namespace

bool evaluate_cylindrical_bell_modes(
    const CylindricalBellModel& model,
    const ModalInfluence& influence,
    const CylindricalModalEvaluationConfig& config,
    ModalData& result) {
    result.clear();
    if (model.mesh_si.nodes.empty() || model.mesh_si.elements.empty()) return false;
    if (influence.radial.size() != model.mesh_si.nodes.size() ||
        influence.axial.size() != model.mesh_si.nodes.size()) return false;
    if (config.harmonics.empty() || config.harmonics.front() != 0 ||
        config.requested_modes <= 0 || config.requested_modes > std::numeric_limits<int>::max() - 8 ||
        config.max_iterations <= 0 ||
        !(config.tolerance > 0.0) || !std::isfinite(config.tolerance) ||
        !(config.minimum_frequency_hz >= 0.0) || !std::isfinite(config.minimum_frequency_hz)) return false;

    const auto system = assemble_cylindrical_bell_system(
        model, config.harmonics, config.angular_samples);
    const auto axis_transform = make_cylindrical_axis_regularity_transform(
        model.mesh_si, system.dofs_per_node, config.harmonics);
    // The free shell retains six rigid-body modes under exact axis regularity.
    // Request one additional complete modal pair as well, so filtering the
    // rigid modes or a frequency-cutoff boundary cannot split a circular pair.
    const int solve_mode_count = config.requested_modes + 8;
    const auto solved = solve_cylindrical_global_modes(
        system, axis_transform, solve_mode_count, config.max_iterations,
        config.tolerance, config.use_rcm_ordering);
    if (solved.frequencies_hz.empty() || !std::isfinite(solved.relative_residual) ||
        solved.relative_residual > std::max(config.tolerance * 10.0, 1.0e-5)) return false;

    // Build the n=0 meridional generalized influence. Higher Fourier basis
    // coefficients are intentionally zero because the current excitation
    // contract is not yet a Chinese-bell strike model.
    const std::size_t basis_count = system.basis_count;
    const std::size_t dpn = system.dofs_per_node;
    std::vector<double> q(system.stiffness.n, 0.0);
    for (std::size_t node = 0; node < model.mesh_si.nodes.size(); ++node) {
        const std::size_t base = node * dpn;
        q[base + 0u * basis_count + 0u] = influence.radial[node];
        q[base + 2u * basis_count + 0u] = influence.axial[node];
    }

    double max_excitation = 0.0;
    std::vector<double> raw_amplitudes;
    raw_amplitudes.reserve(solved.mode_shapes.size());

    for (std::size_t i = 0; i < solved.mode_shapes.size(); ++i) {
        const auto& phi = solved.mode_shapes[i];
        const auto mphi = system.mass.multiply(phi);
        double modal_mass = 0.0;
        double projection = 0.0;
        for (std::size_t j = 0; j < phi.size(); ++j) {
            modal_mass += phi[j] * mphi[j];
            projection += q[j] * mphi[j];
        }
        if (!finite_positive(modal_mass) || solved.frequencies_hz[i] < config.minimum_frequency_hz) continue;
        const double norm = 1.0 / std::sqrt(modal_mass);
        const double excitation = std::abs(projection * norm);
        max_excitation = std::max(max_excitation, excitation);
        raw_amplitudes.push_back(excitation);

        ModalMode mode;
        mode.mode_number = static_cast<int>(result.modes.size()) + 2;
        mode.eigenvalue = solved.eigenvalues[i];
        mode.frequency_hz = solved.frequencies_hz[i];
        mode.modal_mass = 1.0;
        mode.eigen_residual_relative = solved.relative_residual;
        mode.modal_excitation = excitation;
        mode.effective_mass_norm = excitation * excitation;
        mode.amplitude = excitation;
        mode.part_rz = projection * norm;
        mode.mass_rz = mode.part_rz * mode.part_rz;
        mode.fem_displacement.resize(phi.size());
        mode.fem_mass_lumped_weights.assign(phi.size(), 0.0);
        mode.fem_node_r.reserve(model.mesh_si.nodes.size());
        mode.fem_node_z.reserve(model.mesh_si.nodes.size());
        mode.fem_harmonics.clear();
        for (int harmonic : config.harmonics) {
            mode.fem_harmonics.push_back(harmonic);
            if (harmonic != 0) mode.fem_harmonics.push_back(-harmonic);
        }
        for (std::size_t node = 0; node < model.mesh_si.nodes.size(); ++node) {
            mode.fem_node_r.push_back(model.mesh_si.nodes[node].r);
            mode.fem_node_z.push_back(model.mesh_si.nodes[node].z);
        }
        // Row-sum lumping preserves total mass for a consistent mass matrix
        // with nonnegative scalar shape functions. Keep all DOF weights so the
        // tracker can compare complete vector fields in a mass-weighted norm.
        for (int row = 0; row < system.mass.n; ++row) {
            double row_sum = 0.0;
            for (std::size_t k = system.mass.row_offsets[static_cast<std::size_t>(row)];
                 k < system.mass.row_offsets[static_cast<std::size_t>(row) + 1u]; ++k)
                row_sum += system.mass.values[k];
            mode.fem_mass_lumped_weights[static_cast<std::size_t>(row)] = row_sum;
            mode.fem_displacement[static_cast<std::size_t>(row)] = phi[static_cast<std::size_t>(row)] * norm;
        }
        mode.nodal_radial.resize(model.mesh_si.nodes.size(), 0.0);
        mode.nodal_axial.resize(model.mesh_si.nodes.size(), 0.0);
        for (std::size_t node = 0; node < model.mesh_si.nodes.size(); ++node) {
            const std::size_t base = node * dpn;
            mode.nodal_radial[node] = phi[base + 0u * basis_count + 0u] * norm;
            mode.nodal_axial[node] = phi[base + 2u * basis_count + 0u] * norm;
        }
        // The production BellGenerator mesh stores the outer profile first,
        // followed by the inner profile. Export the actual outer-wall samples
        // so shape-based modal continuation works on direct vector-FEM data
        // (and can resample them when the profile mesh is refined).
        const std::size_t outer_count = std::min(model.mesh_si.outer_profile.size(),
                                                 model.mesh_si.nodes.size());
        mode.outer_surface_radial.resize(outer_count);
        mode.outer_surface_axial.resize(outer_count);
        for (std::size_t k = 0; k < outer_count; ++k) {
            mode.outer_surface_radial[k] = mode.nodal_radial[k];
            mode.outer_surface_axial[k] = mode.nodal_axial[k];
        }
        mode.outer_surface_harmonics.clear();
        for (int harmonic : config.harmonics) {
            mode.outer_surface_harmonics.push_back(harmonic); // cosine component
            if (harmonic != 0) mode.outer_surface_harmonics.push_back(-harmonic); // sine component
        }
        mode.outer_surface_fourier_radial.resize(outer_count * basis_count, 0.0);
        mode.outer_surface_fourier_azimuthal.resize(outer_count * basis_count, 0.0);
        mode.outer_surface_fourier_axial.resize(outer_count * basis_count, 0.0);
        for (std::size_t k = 0; k < outer_count; ++k) {
            for (std::size_t b = 0; b < basis_count; ++b) {
                const std::size_t node_base = k * dpn;
                const std::size_t dst = k * basis_count + b;
                mode.outer_surface_fourier_radial[dst] =
                    phi[node_base + 0u * basis_count + b] * norm;
                mode.outer_surface_fourier_azimuthal[dst] =
                    phi[node_base + 1u * basis_count + b] * norm;
                mode.outer_surface_fourier_axial[dst] =
                    phi[node_base + 2u * basis_count + b] * norm;
            }
        }

        // The production mesh stores the inner wall immediately after the outer
        // wall. Export the true inner boundary DOFs rather than estimating them
        // from thickness or offsetting the outer modal field.
        const std::size_t inner_offset = outer_count;
        const std::size_t available_inner = model.mesh_si.nodes.size() > inner_offset
            ? model.mesh_si.nodes.size() - inner_offset : 0u;
        const std::size_t inner_count = std::min(model.mesh_si.inner_profile.size(), available_inner);
        mode.inner_surface_radial.resize(inner_count);
        mode.inner_surface_axial.resize(inner_count);
        mode.inner_surface_harmonics = mode.outer_surface_harmonics;
        mode.inner_surface_fourier_radial.resize(inner_count * basis_count, 0.0);
        mode.inner_surface_fourier_azimuthal.resize(inner_count * basis_count, 0.0);
        mode.inner_surface_fourier_axial.resize(inner_count * basis_count, 0.0);
        for (std::size_t k = 0; k < inner_count; ++k) {
            const std::size_t node_base = (inner_offset + k) * dpn;
            mode.inner_surface_radial[k] = phi[node_base + 0u * basis_count] * norm;
            mode.inner_surface_axial[k] = phi[node_base + 2u * basis_count] * norm;
            for (std::size_t b = 0; b < basis_count; ++b) {
                const std::size_t dst = k * basis_count + b;
                mode.inner_surface_fourier_radial[dst] =
                    phi[node_base + 0u * basis_count + b] * norm;
                mode.inner_surface_fourier_azimuthal[dst] =
                    phi[node_base + 1u * basis_count + b] * norm;
                mode.inner_surface_fourier_axial[dst] =
                    phi[node_base + 2u * basis_count + b] * norm;
            }
        }
        result.modes.push_back(std::move(mode));
    }

    if (result.modes.empty()) return false;
    if (max_excitation > 1.0e-14) {
        for (auto& mode : result.modes) {
            mode.amplitude /= max_excitation;
            mode.modal_excitation /= max_excitation;
        }
    } else {
        for (auto& mode : result.modes) mode.amplitude = 0.0;
    }
    result.fundamental_frequency = result.modes.front().frequency_hz;
    result.sample_rate = 48000.0;
    result.duration = 3.0;
    result.t60_base = 2.0;
    result.target_rms = 0.70;
    result.material = "direct_cylindrical_fem";
    return result.modes.size() >= 5;
}

} // namespace BellOptimizer::FEM

// Chinese-bell strike diagnostics are deliberately kept below the existing
// evaluator. They share the same K/M solve but use a position-dependent
// Fourier load influence rather than the legacy meridional proxy.
namespace BellOptimizer::FEM {


} // namespace BellOptimizer::FEM
