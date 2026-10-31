#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "core/bell_types.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace BellOptimizer;
    BellGeometry g;
    g.H = 1000.0;
    g.R_mouth = 500.0;
    g.t_top = 10.0;
    g.t_bottom = 30.0;
    g.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    g.ellipticity = 0.85;
    BellMaterial mat;
    BellMeshConfig mesh;
    mesh.nNodes = 15;

    const auto model = FEM::make_cylindrical_bell_model(g, mat, mesh);
    FEM::CylindricalModalEvaluationConfig cfg;
    cfg.harmonics = {0, 2, 4};
    cfg.angular_samples = 64;
    cfg.requested_modes = 16;
    cfg.minimum_frequency_hz = 1.0;

    FEM::ChineseBellStrikePositionResult r;
    if (!FEM::optimize_chinese_bell_strike_pose(model, 0.045, 0.22, 5, 12, cfg, r))
        throw std::runtime_error("strike position search failed");
    if (!r.valid || !(r.best_long_selectivity > 0.0) || !(r.best_short_selectivity > 0.0))
        throw std::runtime_error("invalid strike position result");
    if (!std::isfinite(r.best_long_z) || !std::isfinite(r.best_short_z) ||
        !std::isfinite(r.best_long_azimuth_rad) || !std::isfinite(r.best_short_azimuth_rad))
        throw std::runtime_error("non-finite strike position");
    if (!(r.long_target_frequency_hz > 0.0) || !(r.short_target_frequency_hz > 0.0))
        throw std::runtime_error("invalid target frequencies");
    if (!r.family_tracking_valid || !(r.long_family_similarity > 0.60) || !(r.short_family_similarity > 0.60))
        throw std::runtime_error("tracked Chinese modal family is not stable");

    std::cout << "chinese_bell_strike_position_search: PASS\n";
    std::cout << "long_target_mode=" << r.long_target_mode
              << " short_target_mode=" << r.short_target_mode << "\n";
    std::cout << "long_target_frequency_hz=" << r.long_target_frequency_hz
              << " short_target_frequency_hz=" << r.short_target_frequency_hz << "\n";
    std::cout << "family_similarity_long=" << r.long_family_similarity
              << " family_similarity_short=" << r.short_family_similarity
              << " canonical_branch_separated=" << (r.canonical_branch_separated ? 1 : 0)
              << " two_tone_excitation_gate_passed=" << (r.two_tone_excitation_gate_passed ? 1 : 0) << "\n";
    std::cout << "axis_selectivity_long=" << r.long_axis_selectivity
              << " best_long=" << r.best_long_selectivity
              << " z=" << r.best_long_z
              << " phi=" << r.best_long_azimuth_rad << "\n";
    std::cout << "axis_selectivity_short=" << r.short_axis_selectivity
              << " best_short=" << r.best_short_selectivity
              << " z=" << r.best_short_z
              << " phi=" << r.best_short_azimuth_rad << "\n";
    return 0;
}
