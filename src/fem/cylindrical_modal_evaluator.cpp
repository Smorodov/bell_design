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

std::vector<int> production_axis_constraints(const BellMeshData& mesh,
                                             std::size_t dofs_per_node,
                                             std::size_t basis_count) {
    std::vector<int> fixed;
    // The current production bridge uses even harmonics for an elliptic section.
    // At the symmetry-axis node, n=0 radial displacement is zero; higher Fourier
    // displacement amplitudes are pinned at the coordinate singularity. This is
    // a numerical regularity condition, not a physical clamp of the bell wall.
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        if (std::abs(mesh.nodes[node].r) >= 1.0e-9) continue;
        fixed.push_back(static_cast<int>(node * dofs_per_node));
        for (std::size_t b = 1; b < basis_count; ++b) {
            const std::size_t base = node * dofs_per_node + 3u * b;
            fixed.push_back(static_cast<int>(base));
            fixed.push_back(static_cast<int>(base + 1u));
            fixed.push_back(static_cast<int>(base + 2u));
        }
    }
    return fixed;
}

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
        config.requested_modes <= 0 || config.max_iterations <= 0 ||
        !(config.tolerance > 0.0) || !std::isfinite(config.tolerance) ||
        !(config.minimum_frequency_hz >= 0.0) || !std::isfinite(config.minimum_frequency_hz)) return false;

    const auto system = assemble_cylindrical_bell_system(
        model, config.harmonics, config.angular_samples);
    const auto fixed = production_axis_constraints(
        model.mesh_si, system.dofs_per_node, system.basis_count);
    const auto solved = solve_cylindrical_global_modes(
        system, fixed, config.requested_modes, config.max_iterations,
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
        q[base + 0u] = influence.radial[node];
        q[base + 2u] = influence.axial[node];
    }
    (void)basis_count;

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
            mode.nodal_radial[node] = phi[base + 0u] * norm;
            mode.nodal_axial[node] = phi[base + 2u] * norm;
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
                const std::size_t dof = k * dpn + 3u * b;
                const std::size_t dst = k * basis_count + b;
                mode.outer_surface_fourier_radial[dst] = phi[dof] * norm;
                mode.outer_surface_fourier_azimuthal[dst] = phi[dof + 1u] * norm;
                mode.outer_surface_fourier_axial[dst] = phi[dof + 2u] * norm;
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

bool evaluate_chinese_bell_strike_selectivity(
    const CylindricalBellModel& model,
    double center_z,
    double meridional_width,
    double angular_width_rad,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikeSelectivityResult& result) {
    result = {};
    if (model.mesh_si.nodes.empty() || model.mesh_si.elements.empty()) return false;
    if (!std::isfinite(center_z) || !(meridional_width > 0.0) || !(angular_width_rad > 0.0)) return false;

    const auto system = assemble_cylindrical_bell_system(
        model, config.harmonics, config.angular_samples);
    const auto fixed = production_axis_constraints(
        model.mesh_si, system.dofs_per_node, system.basis_count);
    const auto solved = solve_cylindrical_global_modes(
        system, fixed, config.requested_modes, config.max_iterations,
        config.tolerance, config.use_rcm_ordering);
    if (solved.frequencies_hz.empty() || !std::isfinite(solved.relative_residual) ||
        solved.relative_residual > std::max(config.tolerance * 10.0, 1.0e-5)) return false;

    Excitation::ChineseBellStrikeConfig long_cfg;
    long_cfg.center_z = center_z;
    long_cfg.azimuth_rad = 0.0;
    long_cfg.meridional_width = meridional_width;
    long_cfg.angular_width_rad = angular_width_rad;
    Excitation::ChineseBellStrikeConfig short_cfg = long_cfg;
    short_cfg.azimuth_rad = 0.5 * PI;

    const auto long_influence = Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, config.harmonics, model.ellipticity, long_cfg, config.angular_samples);
    const auto short_influence = Excitation::make_chinese_bell_strike_influence(
        model.mesh_si, config.harmonics, model.ellipticity, short_cfg, config.angular_samples);
    if (long_influence.dof_values.size() != static_cast<std::size_t>(system.stiffness.n) ||
        short_influence.dof_values.size() != static_cast<std::size_t>(system.stiffness.n)) return false;

    result.frequencies_hz.reserve(solved.frequencies_hz.size());
    result.long_axis_excitation.reserve(solved.mode_shapes.size());
    result.short_axis_excitation.reserve(solved.mode_shapes.size());

    double max_long = 0.0, max_short = 0.0;
    for (std::size_t i = 0; i < solved.mode_shapes.size(); ++i) {
        if (solved.frequencies_hz[i] < config.minimum_frequency_hz) continue;
        const auto& phi = solved.mode_shapes[i];
        const auto mphi = system.mass.multiply(phi);
        double modal_mass = 0.0, ql = 0.0, qs = 0.0;
        for (std::size_t j = 0; j < phi.size(); ++j) {
            modal_mass += phi[j] * mphi[j];
            // FEMStrikeInfluence is a nodal generalized-force distribution.
            // Modal force is phi^T f; multiplying f by M here would turn the
            // force into a mass-weighted proxy and make the result mesh/material
            // dependent in the wrong way.
            ql += long_influence.dof_values[j] * phi[j];
            qs += short_influence.dof_values[j] * phi[j];
        }
        if (!(modal_mass > 0.0) || !std::isfinite(modal_mass)) continue;
        const double scale = 1.0 / std::sqrt(modal_mass);
        const double el = std::abs(ql * scale);
        const double es = std::abs(qs * scale);
        result.frequencies_hz.push_back(solved.frequencies_hz[i]);
        result.long_axis_excitation.push_back(el);
        result.short_axis_excitation.push_back(es);
        max_long = std::max(max_long, el);
        max_short = std::max(max_short, es);
    }
    if (result.frequencies_hz.size() < 2 || !(max_long > 0.0) || !(max_short > 0.0)) return false;
    for (double& v : result.long_axis_excitation) v /= max_long;
    for (double& v : result.short_axis_excitation) v /= max_short;

    // For each canonical strike, measure how much of the strongest competing
    // response remains outside the dominant mode. This is a diagnostic only;
    // later Chinese-bell family tracking will replace rank-based selection.
    const auto max_index = [](const std::vector<double>& v) {
        return static_cast<std::size_t>(std::distance(v.begin(), std::max_element(v.begin(), v.end())));
    };
    const std::size_t il = max_index(result.long_axis_excitation);
    const std::size_t is = max_index(result.short_axis_excitation);
    double second_l = 0.0, second_s = 0.0;
    for (std::size_t i = 0; i < result.frequencies_hz.size(); ++i) {
        if (i != il) second_l = std::max(second_l, result.long_axis_excitation[i]);
        if (i != is) second_s = std::max(second_s, result.short_axis_excitation[i]);
    }
    result.long_axis_cross_ratio = second_l;
    result.short_axis_cross_ratio = second_s;
    result.valid = true;
    return true;
}


bool optimize_chinese_bell_strike_pose(
    const CylindricalBellModel& model,
    double meridional_width,
    double angular_width_rad,
    unsigned z_samples,
    unsigned azimuth_samples,
    const CylindricalModalEvaluationConfig& config,
    ChineseBellStrikePositionResult& result) {
    result = {};
    if (model.mesh_si.nodes.empty() || model.mesh_si.elements.empty() ||
        !(meridional_width > 0.0) || !(angular_width_rad > 0.0) ||
        !std::isfinite(meridional_width) || !std::isfinite(angular_width_rad) ||
        z_samples < 2 || azimuth_samples < 4) return false;

    const auto system = assemble_cylindrical_bell_system(
        model, config.harmonics, config.angular_samples);
    const auto fixed = production_axis_constraints(
        model.mesh_si, system.dofs_per_node, system.basis_count);
    const auto solved = solve_cylindrical_global_modes(
        system, fixed, config.requested_modes, config.max_iterations,
        config.tolerance, config.use_rcm_ordering);
    if (solved.frequencies_hz.empty() || !std::isfinite(solved.relative_residual) ||
        solved.relative_residual > std::max(config.tolerance * 10.0, 1.0e-5)) return false;

    // Build a numerically near-circular reference at the same mesh and material.
    // A tiny ellipticity regularizes the exact repeated eigenspace in the current
    // sparse solver; the first non-axisymmetric pair is tracked as a subspace;
    // at q=1 its individual eigenvectors are basis-dependent, so the family
    // is tracked as a subspace rather than by sorted frequency rank.
    auto circular_reference_model = model;
    circular_reference_model.ellipticity = 0.9999;
    const auto reference_system = assemble_cylindrical_bell_system(
        circular_reference_model, config.harmonics, config.angular_samples);
    const auto reference_fixed = production_axis_constraints(
        model.mesh_si, reference_system.dofs_per_node, reference_system.basis_count);
    const auto reference = solve_cylindrical_global_modes(
        reference_system, reference_fixed, std::max(config.requested_modes, 8),
        std::max(config.max_iterations * 2, 2400), config.tolerance, config.use_rcm_ordering);
    if (reference.mode_shapes.size() < 2 ||
        !std::isfinite(reference.relative_residual) ||
        reference.relative_residual > std::max(config.tolerance * 10.0, 1.0e-5)) return false;

    // Seed the family with the lowest-frequency nearly-degenerate pair in the
    // circular reference.  Sorted rank alone is not a family identity, so the
    // search explicitly verifies degeneracy before using the pair as a tracking
    // subspace.
    std::size_t ref_i = 0, ref_j = 1;
    double best_ref_split = std::numeric_limits<double>::infinity();
    const std::size_t ref_limit = std::min<std::size_t>(reference.frequencies_hz.size(), 12);
    for (std::size_t i = 0; i < ref_limit; ++i) {
        for (std::size_t j = i + 1; j < ref_limit; ++j) {
            const double split = std::abs(reference.frequencies_hz[j] - reference.frequencies_hz[i]) /
                std::max(reference.frequencies_hz[i], 1.0e-300);
            if (split < best_ref_split) { best_ref_split = split; ref_i = i; ref_j = j; }
        }
    }
    if (!(best_ref_split < 1.0e-3)) return false;

    const auto& family_metric = reference_system.mass;
    ChineseBellStrikeSelectivityResult canonical;
    if (!evaluate_chinese_bell_strike_selectivity(
            model, 0.70 * (model.mesh_si.nodes.back().z), meridional_width,
            angular_width_rad, config, canonical) || !canonical.valid) return false;
    if (canonical.frequencies_hz.size() < 2) return false;
    const auto max_index = [](const std::vector<double>& v) {
        return static_cast<std::size_t>(std::distance(v.begin(),
            std::max_element(v.begin(), v.end())));
    };
    const std::size_t long_target = max_index(canonical.long_axis_excitation);
    const std::size_t short_target = max_index(canonical.short_axis_excitation);
    if (long_target >= solved.mode_shapes.size() || short_target >= solved.mode_shapes.size()) return false;

    std::vector<std::vector<double>> normalized_shapes;
    normalized_shapes.reserve(solved.mode_shapes.size());
    for (const auto& phi : solved.mode_shapes) {
        const auto mphi = system.mass.multiply(phi);
        double mm = 0.0;
        for (std::size_t j = 0; j < phi.size(); ++j) mm += phi[j] * mphi[j];
        if (!(mm > 0.0) || !std::isfinite(mm)) return false;
        const double inv = 1.0 / std::sqrt(mm);
        std::vector<double> nphi(phi.size());
        // Unit-modal-mass displacement shape. Generalized force is evaluated
        // later as f^T phi/sqrt(phi^T M phi), not f^T M phi.
        for (std::size_t j = 0; j < phi.size(); ++j) nphi[j] = phi[j] * inv;
        normalized_shapes.push_back(std::move(nphi));
    }

    const double z_min = model.mesh_si.outer_profile.front().z;
    const double z_max = model.mesh_si.outer_profile.back().z;
    if (!std::isfinite(z_min) || !std::isfinite(z_max) || !(z_max > z_min)) return false;

    auto selectivity = [&](double z, double azimuth, std::size_t target) {
        Excitation::ChineseBellStrikeConfig strike;
        strike.center_z = z;
        strike.azimuth_rad = azimuth;
        strike.meridional_width = meridional_width;
        strike.angular_width_rad = angular_width_rad;
        const auto influence = Excitation::make_chinese_bell_strike_influence(
            model.mesh_si, config.harmonics, model.ellipticity, strike, config.angular_samples);
        if (influence.dof_values.size() != normalized_shapes[target].size()) return 0.0;
        double target_q = 0.0;
        double competitor_sq = 0.0;
        for (std::size_t i = 0; i < normalized_shapes.size(); ++i) {
            double q = 0.0;
            for (std::size_t j = 0; j < influence.dof_values.size(); ++j)
                q += influence.dof_values[j] * normalized_shapes[i][j];
            q = std::abs(q);
            if (i == target) target_q = q;
            else competitor_sq += q * q;
        }
        const double denom = std::sqrt(competitor_sq) + 1.0e-15;
        return target_q / denom;
    };

    result.long_target_mode = long_target;
    result.short_target_mode = short_target;
    result.long_target_frequency_hz = solved.frequencies_hz[long_target];
    result.short_target_frequency_hz = solved.frequencies_hz[short_target];
    result.long_axis_selectivity = selectivity(z_min + 0.70 * (z_max - z_min), 0.0, long_target);
    result.short_axis_selectivity = selectivity(z_min + 0.70 * (z_max - z_min), 0.5 * PI, short_target);

    // Track the physical two-dimensional family into the elliptic candidate.
    // The family metric is deliberately the circular reference mass matrix: it
    // is a fixed tracking coordinate, not an extra q-dependent objective term.
    double best_family_similarity = 0.0;
    std::size_t family_i = 0, family_j = 1;
    for (std::size_t i = 0; i < solved.mode_shapes.size(); ++i) {
        for (std::size_t j = i + 1; j < solved.mode_shapes.size(); ++j) {
            const double sim = mass_weighted_subspace_similarity(
                reference, {ref_i, ref_j}, solved, {i, j}, family_metric);
            if (sim > best_family_similarity) {
                best_family_similarity = sim;
                family_i = i;
                family_j = j;
            }
        }
    }
    result.long_family_similarity = best_family_similarity;
    result.short_family_similarity = best_family_similarity;
    result.family_tracking_valid = best_family_similarity > 0.60;
    if (!result.family_tracking_valid) return false;

    // The tracked family, not an arbitrary sorted mode, is now the provisional
    // Chinese two-tone subspace.  The two canonical strikes must select opposite
    // split branches; this is the first explicit cross-excitation contract.
    auto branch_projections = [&](double azimuth) {
        Excitation::ChineseBellStrikeConfig strike;
        strike.center_z = z_min + 0.70 * (z_max - z_min);
        strike.azimuth_rad = azimuth;
        strike.meridional_width = meridional_width;
        strike.angular_width_rad = angular_width_rad;
        const auto influence = Excitation::make_chinese_bell_strike_influence(
            model.mesh_si, config.harmonics, model.ellipticity, strike, config.angular_samples);
        auto projection = [&](std::size_t mode) {
            const auto& phi = solved.mode_shapes[mode];
            const auto mphi = system.mass.multiply(phi);
            double mm = 0.0, q = 0.0;
            for (std::size_t k = 0; k < phi.size(); ++k) {
                mm += phi[k] * mphi[k];
                q += influence.dof_values[k] * phi[k];
            }
            return mm > 0.0 ? std::abs(q) / std::sqrt(mm) : 0.0;
        };
        return std::pair<double, double>{projection(family_i), projection(family_j)};
    };
    const auto long_branches = branch_projections(0.0);
    const auto short_branches = branch_projections(0.5 * PI);
    const double long_contrast = long_branches.first / (long_branches.second + 1.0e-15);
    const double long_contrast_swap = long_branches.second / (long_branches.first + 1.0e-15);
    const double short_contrast = short_branches.first / (short_branches.second + 1.0e-15);
    const double short_contrast_swap = short_branches.second / (short_branches.first + 1.0e-15);
    const bool orientation_a = long_contrast >= long_contrast_swap;
    const bool orientation_b = short_contrast >= short_contrast_swap;
    // Do not manufacture two-tone behavior when both canonical strikes prefer
    // the same branch.  Keep the family tracking result valid, but expose the
    // failure as a separate diagnostic so the next objective stage can reject
    // the geometry for insufficient excitation separation.
    result.canonical_branch_separated = (orientation_a != orientation_b);
    result.two_tone_excitation_gate_passed = result.canonical_branch_separated;
    const std::size_t tracked_long = orientation_a ? family_i : family_j;
    const std::size_t tracked_short = orientation_a ? family_j : family_i;
    result.long_target_mode = tracked_long;
    result.short_target_mode = tracked_short;
    result.long_target_frequency_hz = solved.frequencies_hz[tracked_long];
    result.short_target_frequency_hz = solved.frequencies_hz[tracked_short];

    result.best_long_selectivity = -1.0;
    result.best_short_selectivity = -1.0;

    for (unsigned iz = 0; iz < z_samples; ++iz) {
        const double z = z_min + (z_max - z_min) * static_cast<double>(iz) /
            static_cast<double>(z_samples - 1);
        for (unsigned ia = 0; ia < azimuth_samples; ++ia) {
            const double azimuth = PI * static_cast<double>(ia) /
                static_cast<double>(azimuth_samples);
            const double sl = selectivity(z, azimuth, tracked_long);
            const double ss = selectivity(z, azimuth, tracked_short);
            if (std::isfinite(sl) && sl > result.best_long_selectivity) {
                result.best_long_selectivity = sl;
                result.best_long_z = z;
                result.best_long_azimuth_rad = azimuth;
            }
            if (std::isfinite(ss) && ss > result.best_short_selectivity) {
                result.best_short_selectivity = ss;
                result.best_short_z = z;
                result.best_short_azimuth_rad = azimuth;
            }
        }
    }
    result.valid = result.best_long_selectivity >= 0.0 && result.best_short_selectivity >= 0.0;
    return result.valid;
}

} // namespace BellOptimizer::FEM
