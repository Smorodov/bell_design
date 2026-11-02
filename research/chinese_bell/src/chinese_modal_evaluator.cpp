#include "chinese_modal_evaluator.hpp"
#include "excitation/chinese_bell_strike.hpp"
#include "fem/cylindrical_global_modes.hpp"
#include "fem/cylindrical_modal_tracking.hpp"
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

} // namespace

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
    result.mode_indices.reserve(solved.frequencies_hz.size());
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
        result.mode_indices.push_back(i);
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
    ChineseBellStrikePositionResult& result,
    double canonical_strike_z_fraction) {
    result = {};
    if (model.mesh_si.nodes.empty() || model.mesh_si.elements.empty() ||
        !(meridional_width > 0.0) || !(angular_width_rad > 0.0) ||
        !std::isfinite(meridional_width) || !std::isfinite(angular_width_rad) ||
        !std::isfinite(canonical_strike_z_fraction) || canonical_strike_z_fraction < 0.0 ||
        canonical_strike_z_fraction > 1.0 || z_samples < 2 || azimuth_samples < 4) return false;

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
            model,
            model.mesh_si.outer_profile.front().z + canonical_strike_z_fraction *
                (model.mesh_si.outer_profile.back().z - model.mesh_si.outer_profile.front().z),
            meridional_width,
            angular_width_rad, config, canonical) || !canonical.valid) return false;
    if (canonical.frequencies_hz.size() < 2) return false;
    const auto max_index = [](const std::vector<double>& v) {
        return static_cast<std::size_t>(std::distance(v.begin(),
            std::max_element(v.begin(), v.end())));
    };
    const std::size_t long_target_position = max_index(canonical.long_axis_excitation);
    const std::size_t short_target_position = max_index(canonical.short_axis_excitation);
    if (long_target_position >= canonical.mode_indices.size() ||
        short_target_position >= canonical.mode_indices.size()) return false;
    const std::size_t long_target = canonical.mode_indices[long_target_position];
    const std::size_t short_target = canonical.mode_indices[short_target_position];
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
    const double canonical_z = z_min + canonical_strike_z_fraction * (z_max - z_min);
    result.long_axis_selectivity = selectivity(canonical_z, 0.0, long_target);
    result.short_axis_selectivity = selectivity(canonical_z, 0.5 * PI, short_target);

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
        strike.center_z = z_min + canonical_strike_z_fraction * (z_max - z_min);
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
    result.canonical_branch_orientation = orientation_a == orientation_b ? 0 : (orientation_a ? 1 : -1);
    result.long_branch_contrast = std::max(long_contrast, long_contrast_swap);
    result.short_branch_contrast = std::max(short_contrast, short_contrast_swap);
    result.two_tone_excitation_gate_passed = result.canonical_branch_separated;
    const std::size_t tracked_long = orientation_a ? family_i : family_j;
    const std::size_t tracked_short = orientation_a ? family_j : family_i;
    result.long_target_mode = tracked_long;
    result.short_target_mode = tracked_short;
    result.long_target_frequency_hz = solved.frequencies_hz[tracked_long];
    result.short_target_frequency_hz = solved.frequencies_hz[tracked_short];
    result.pair_frequency_split_cents = 1200.0 * std::log2(
        result.short_target_frequency_hz / result.long_target_frequency_hz);

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

