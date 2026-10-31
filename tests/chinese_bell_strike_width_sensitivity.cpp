#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "core/bell_types.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using namespace BellOptimizer;
    BellGeometry geometry;
    geometry.H = 1000.0;
    geometry.R_mouth = 500.0;
    geometry.t_top = 10.0;
    geometry.t_bottom = 30.0;
    geometry.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    geometry.ellipticity = 0.85;
    BellMaterial material;
    BellMeshConfig mesh;
    mesh.nNodes = 31;
    const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh);

    FEM::CylindricalModalEvaluationConfig config;
    config.harmonics = {0, 2, 4};
    config.angular_samples = 64;
    config.requested_modes = 16;
    config.minimum_frequency_hz = 1.0;

    const std::vector<double> widths{0.015, 0.035, 0.070};
    std::vector<FEM::ChineseBellStrikeSelectivityResult> results;
    for (double width : widths) {
        FEM::ChineseBellStrikeSelectivityResult result;
        if (!FEM::evaluate_chinese_bell_strike_selectivity(
                model, 0.62, width, 0.20, config, result))
            throw std::runtime_error("strike width evaluation failed");
        if (result.frequencies_hz.size() < 5 ||
            result.long_axis_excitation.size() != result.frequencies_hz.size() ||
            result.short_axis_excitation.size() != result.frequencies_hz.size())
            throw std::runtime_error("inconsistent strike-width result dimensions");
        for (std::size_t i = 0; i < result.frequencies_hz.size(); ++i) {
            if (!(result.frequencies_hz[i] > 0.0) ||
                !std::isfinite(result.long_axis_excitation[i]) ||
                !std::isfinite(result.short_axis_excitation[i]))
                throw std::runtime_error("non-finite modal participation");
        }
        results.push_back(std::move(result));
    }

    const auto& reference = results.front();
    double max_frequency_delta = 0.0;
    double max_width_response_delta = 0.0;
    for (std::size_t w = 1; w < results.size(); ++w) {
        if (results[w].frequencies_hz.size() != reference.frequencies_hz.size())
            throw std::runtime_error("mode count changed with contact width");
        for (std::size_t i = 0; i < reference.frequencies_hz.size(); ++i) {
            max_frequency_delta = std::max(max_frequency_delta,
                std::abs(results[w].frequencies_hz[i] - reference.frequencies_hz[i]));
            max_width_response_delta = std::max(max_width_response_delta,
                std::abs(results[w].long_axis_excitation[i] - reference.long_axis_excitation[i]));
            max_width_response_delta = std::max(max_width_response_delta,
                std::abs(results[w].short_axis_excitation[i] - reference.short_axis_excitation[i]));
        }
    }
    if (max_frequency_delta > 1.0e-8)
        throw std::runtime_error("contact width changed the structural eigenfrequencies");
    if (!(max_width_response_delta > 1.0e-5))
        throw std::runtime_error("modal participation is insensitive to contact width");

    std::cout << "chinese_bell_strike_width_sensitivity: PASS\n"
              << "max_frequency_delta_hz=" << max_frequency_delta << '\n'
              << "max_normalized_participation_delta=" << max_width_response_delta << '\n';
    return 0;
}
