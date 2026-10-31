#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "core/bell_types.h"
#include <algorithm>
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
    cfg.requested_modes = 20;
    cfg.minimum_frequency_hz = 1.0;

    FEM::ChineseBellStrikeSelectivityResult r;
    const double center_z = 0.70;
    if (!FEM::evaluate_chinese_bell_strike_selectivity(model, center_z, 0.045, 0.22, cfg, r))
        throw std::runtime_error("Chinese strike selectivity evaluation failed");
    if (!r.valid || r.frequencies_hz.size() < 5) throw std::runtime_error("insufficient strike modes");
    if (r.long_axis_excitation.size() != r.frequencies_hz.size() ||
        r.short_axis_excitation.size() != r.frequencies_hz.size())
        throw std::runtime_error("strike projection size mismatch");
    for (double v : r.long_axis_excitation) if (!std::isfinite(v) || v < 0.0) throw std::runtime_error("invalid long-axis excitation");
    for (double v : r.short_axis_excitation) if (!std::isfinite(v) || v < 0.0) throw std::runtime_error("invalid short-axis excitation");
    const double long_peak = *std::max_element(r.long_axis_excitation.begin(), r.long_axis_excitation.end());
    const double short_peak = *std::max_element(r.short_axis_excitation.begin(), r.short_axis_excitation.end());
    if (std::abs(long_peak - 1.0) > 1e-12 || std::abs(short_peak - 1.0) > 1e-12)
        throw std::runtime_error("strike normalization failed");
    if (!std::isfinite(r.long_axis_cross_ratio) || !std::isfinite(r.short_axis_cross_ratio))
        throw std::runtime_error("invalid cross excitation diagnostic");

    std::cout << "chinese_bell_strike_selectivity: PASS\n";
    std::cout << "modes=" << r.frequencies_hz.size()
              << " long_cross=" << r.long_axis_cross_ratio
              << " short_cross=" << r.short_axis_cross_ratio << "\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(5, r.frequencies_hz.size()); ++i)
        std::cout << r.frequencies_hz[i] << " Hz long=" << r.long_axis_excitation[i]
                  << " short=" << r.short_axis_excitation[i] << "\n";
    return 0;
}
