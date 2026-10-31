#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include <iostream>
#include <vector>

int main() {
    using namespace BellOptimizer;
    BellGeometry g;
    g.H = 1000.0; g.R_mouth = 500.0; g.t_top = 10.0; g.t_bottom = 30.0;
    g.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0}; g.ellipticity = 0.85;
    BellMaterial material;
    struct Case { int nodes; std::vector<int> harmonics; unsigned angular; };
    const std::vector<Case> cases = {
        {15, {0,2,4}, 64}, {23, {0,2,4}, 64},
        {15, {0,2,4,6}, 64}, {15, {0,2,4,6}, 128}
    };
    for (const auto& c : cases) {
        BellMeshConfig mc; mc.nNodes = c.nodes;
        const auto model = FEM::make_cylindrical_bell_model(g, material, mc);
        FEM::CylindricalModalEvaluationConfig cfg;
        cfg.harmonics = c.harmonics; cfg.angular_samples = c.angular;
        cfg.requested_modes = 16; cfg.minimum_frequency_hz = 1.0;
        FEM::ChineseBellStrikePositionResult r;
        if (!FEM::optimize_chinese_bell_strike_pose(model, 0.045, 0.22, 5, 12, cfg, r) || !r.valid) { std::cout << "nodes=" << c.nodes << " harmonics=" << c.harmonics.size() << " angular=" << c.angular << " result=UNAVAILABLE gate=0\n"; continue; }
        std::cout << "nodes=" << c.nodes << " harmonics=" << c.harmonics.size()
                  << " angular=" << c.angular
                  << " f_long=" << r.long_target_frequency_hz
                  << " f_short=" << r.short_target_frequency_hz
                  << " split_cents=" << (1200.0 * std::log2(r.short_target_frequency_hz / r.long_target_frequency_hz))
                  << " family_similarity=" << r.long_family_similarity
                  << " separated=" << (r.canonical_branch_separated ? 1 : 0)
                  << " gate=" << (r.two_tone_excitation_gate_passed ? 1 : 0) << '\n';
    }
    return 0;
}
