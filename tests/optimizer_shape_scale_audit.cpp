#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool solve(const BellGeometry& geometry, BellOptimizer::ModalData& result) {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return false;
    BellOptimizer::InternalModalSolver solver(true, 8);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "shape_scale_audit", &modal_influence_1};
    return solver.solve(problem, result);
}

} // namespace

int main() {
    const BellGeometry base{1000.0, 500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0}};
    BellGeometry scaled = base;
    scaled.scale(1.7);

    BellOptimizer::ModalData a;
    BellOptimizer::ModalData b;
    if (!solve(base, a) || !solve(scaled, b) || a.modes.size() < 5 || b.modes.size() < 5) {
        std::cerr << "modal solve failed\n";
        return 2;
    }

    double max_ratio_error = 0.0;
    double max_amp_error = 0.0;
    for (int i = 0; i < 5; ++i) {
        const double expected_frequency = a.modes[static_cast<size_t>(i)].frequency_hz / 1.7;
        max_ratio_error = std::max(max_ratio_error,
            std::abs((b.modes[static_cast<size_t>(i)].frequency_hz / expected_frequency) - 1.0));
    }
    for (int i = 1; i < 5; ++i) {
        const double ar = a.modes[static_cast<size_t>(i)].amplitude / a.modes[0].amplitude;
        const double br = b.modes[static_cast<size_t>(i)].amplitude / b.modes[0].amplitude;
        max_amp_error = std::max(max_amp_error, std::abs(ar - br));
    }

    std::cout << "max_frequency_scale_error=" << max_ratio_error << '\n';
    std::cout << "max_amplitude_ratio_error=" << max_amp_error << '\n';

    if (!(max_ratio_error < 1e-9)) return 3;
    if (!(max_amp_error < 1e-6)) return 4;
    return 0;
}
