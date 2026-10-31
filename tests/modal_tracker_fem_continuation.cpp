#include "core/modal_state_tracker.hpp"
#include "core/modal_tracking.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static BellGeometry geometry(double q) {
    BellGeometry g;
    g.H = 658.3175709034034;
    g.R_mouth = 355.9758748202059;
    g.t_top = 22.85675347179014;
    g.t_bottom = 46.62524667849419;
    g.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    g.ellipticity = q;
    return g;
}

static ModalData solve(double q, int nodes, const std::vector<int>& harmonics,
                       unsigned angular_samples) {
    BellMaterial material;
    BellMeshConfig mesh;
    mesh.nNodes = nodes;
    auto g = geometry(q);
    g.scale(1000.0 / g.H);
    const auto model = make_cylindrical_bell_model(g, material, mesh);
    CylindricalModalEvaluationConfig cfg;
    cfg.harmonics = harmonics;
    cfg.angular_samples = angular_samples;
    cfg.requested_modes = 14;
    cfg.max_iterations = 5000;
    cfg.tolerance = 1.0e-6;
    cfg.minimum_frequency_hz = 1.0;
    ModalData result;
    Excitation::MeridionalInfluence influence;
    require(evaluate_cylindrical_bell_modes(model, influence.make(model.mesh_si), cfg, result),
            "direct FEM solve failed");
    require(result.modes.size() >= 6, "direct FEM returned too few modes");
    for (const auto& mode : result.modes) {
        require(mode.outer_surface_radial.size() == static_cast<std::size_t>(nodes),
                "direct FEM did not export outer profile samples");
        require(mode.outer_surface_axial.size() == static_cast<std::size_t>(nodes),
                "direct FEM outer profile sample sizes differ");
        std::vector<int> expanded;
        for (int h : harmonics) { expanded.push_back(h); if (h != 0) expanded.push_back(-h); }
        require(mode.outer_surface_harmonics == expanded,
                "direct FEM Fourier metadata was not preserved");
        require(mode.outer_surface_fourier_radial.size() ==
                    static_cast<std::size_t>(nodes) * expanded.size() &&
                mode.outer_surface_fourier_azimuthal.size() ==
                    static_cast<std::size_t>(nodes) * expanded.size() &&
                mode.outer_surface_fourier_axial.size() ==
                    static_cast<std::size_t>(nodes) * expanded.size(),
                "direct FEM Fourier boundary data has invalid size");
        require(mode.fem_node_r.size() == mode.fem_node_z.size() &&
                mode.fem_displacement.size() == mode.fem_node_r.size() * expanded.size() * 3u &&
                mode.fem_mass_lumped_weights.size() == mode.fem_displacement.size() &&
                mode.fem_harmonics == expanded,
                "complete FEM field or lumped-mass metadata has invalid dimensions");
        double weighted_norm = 0.0;
        for (std::size_t k = 0; k < mode.fem_displacement.size(); ++k) {
            require(std::isfinite(mode.fem_displacement[k]) &&
                    std::isfinite(mode.fem_mass_lumped_weights[k]) &&
                    mode.fem_mass_lumped_weights[k] >= -1e-12,
                    "invalid full-field value or mass-lumped weight");
            weighted_norm += std::max(0.0, mode.fem_mass_lumped_weights[k]) *
                             mode.fem_displacement[k] * mode.fem_displacement[k];
        }
        require(weighted_norm > 1e-20, "full FEM field has zero mass-weighted norm");
    }
    return result;
}

int main() {
    try {
        const std::vector<int> base_harmonics{0, 2, 4};
        auto reference = solve(0.9999, 15, base_harmonics, 64);
        ModalStateTracker tracker;
        require(tracker.initialize(reference, 0, 1, 0.9999),
                "Kalman/subspace tracker initialization failed");

        std::size_t kalman_updates = 0;
        std::size_t kalman_ambiguous = 0;
        std::size_t shape_matches = 0;
        double min_shape_similarity = 1.0;
        double min_pair_similarity = 1.0;
        ModalData previous = reference;

        for (double q : {0.99, 0.98, 0.97, 0.96}) {
            auto current = solve(q, 15, base_harmonics, 64);

            // Baseline: global shape assignment, without a frequency-state predictor.
            auto shape_candidate = current;
            const auto shape_report = trackModalFamilies(previous, shape_candidate);
            if (shape_report.used_shape_tracking) {
                ++shape_matches;
                min_shape_similarity = std::min(min_shape_similarity,
                                                shape_report.minimum_similarity);
            }

            // Candidate method: Kalman frequency prediction plus 2-D subspace matching.
            const auto kalman_report = tracker.update(current, q);
            if (kalman_report.updated) {
                ++kalman_updates;
                min_pair_similarity = std::min(min_pair_similarity,
                                               kalman_report.pair_subspace_similarity);
            } else if (kalman_report.ambiguous) {
                ++kalman_ambiguous;
            }

            std::cout << "q=" << q
                      << " shape_match=" << (shape_report.used_shape_tracking ? 1 : 0)
                      << " shape_min=" << shape_report.minimum_similarity
                      << " kalman_update=" << (kalman_report.updated ? 1 : 0)
                      << " kalman_ambiguous=" << (kalman_report.ambiguous ? 1 : 0)
                      << " kalman_pair_similarity=" << kalman_report.pair_subspace_similarity
                      << " kalman_margin=" << kalman_report.assignment_margin
                      << " mass_shape_a=" << kalman_report.tracks[0].shape_similarity
                      << " mass_shape_b=" << kalman_report.tracks[1].shape_similarity << '\n';
            previous = std::move(current);
        }

        // Cross-mesh continuation: 15-node reference to 23-node candidate.
        // This specifically exercises profile-coordinate resampling.
        auto refined = solve(0.95, 23, base_harmonics, 64);
        auto refined_baseline = refined;
        const auto cross_mesh_shape = trackModalFamilies(previous, refined_baseline);
        const auto cross_mesh = tracker.update(refined, 0.95);
        std::cout << "cross_mesh_shape_match=" << (cross_mesh_shape.used_shape_tracking ? 1 : 0)
                  << " cross_mesh update=" << (cross_mesh.updated ? 1 : 0)
                  << " ambiguous=" << (cross_mesh.ambiguous ? 1 : 0)
                  << " pair_similarity=" << cross_mesh.pair_subspace_similarity
                  << " margin=" << cross_mesh.assignment_margin << '\n';

        // Fourier enrichment at fixed mesh: do not assume the state tracker
        // must accept a basis change; report whether the physical family remains
        // identifiable and fail only on malformed output or invalid numerics.
        auto enriched = solve(0.94, 23, {0, 2, 4, 6}, 128);
        auto enriched_baseline = enriched;
        const auto enriched_shape = trackModalFamilies(refined, enriched_baseline);
        const auto enriched_report = tracker.update(enriched, 0.94);
        std::cout << "fourier_enriched_shape_match=" << (enriched_shape.used_shape_tracking ? 1 : 0)
                  << " fourier_enriched update=" << (enriched_report.updated ? 1 : 0)
                  << " ambiguous=" << (enriched_report.ambiguous ? 1 : 0)
                  << " no_compatible=" << (enriched_report.no_compatible_assignment ? 1 : 0)
                  << " pair_similarity=" << enriched_report.pair_subspace_similarity
                  << " margin=" << enriched_report.assignment_margin << '\n';

        require(shape_matches >= 2, "baseline shape matcher did not produce enough FEM matches");
        require(kalman_updates >= 4,
                "state tracker failed to continue the same-mesh FEM sequence");
        require(cross_mesh.updated,
                "state tracker failed to continue across profile mesh refinement");
        require(enriched_report.no_compatible_assignment || enriched_report.ambiguous ||
                    enriched_report.updated,
                "state tracker returned an unclassified Fourier-enrichment result");
        require(std::isfinite(min_shape_similarity) && std::isfinite(min_pair_similarity),
                "non-finite tracker diagnostics");

        std::cout << "summary shape_matches=" << shape_matches
                  << " kalman_updates=" << kalman_updates
                  << " kalman_ambiguous=" << kalman_ambiguous
                  << " min_shape_similarity=" << min_shape_similarity
                  << " min_pair_similarity=" << min_pair_similarity << '\n';
        std::cout << "modal_tracker_fem_continuation: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "modal_tracker_fem_continuation: FAIL: " << e.what() << '\n';
        return 1;
    }
}
