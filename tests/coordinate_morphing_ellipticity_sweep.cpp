#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "math/modal_coordinate_morphing.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

bool solve(const BellGeometry& geometry, const BellMaterial& material,
           BellOptimizer::ModalData& result) {
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return false;
    BellOptimizer::InternalModalSolver solver(true, 6, true);
    BellOptimizer::Excitation::MeridionalInfluence influence;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "ellipticity_sweep", &influence};
    return solver.solve(problem, result);
}

} // namespace

int main() {
    BellMaterial material;
    const BellGeometry geometry{1000.0, 500.0, 10.0, 30.0,
                                {0.0, 0.15, 0.40, 0.70, 1.0}};

    BellOptimizer::ModalData base;
    if (!solve(geometry, material, base) || base.modes.size() < 4) {
        std::cerr << "failed baseline solve\n";
        return 2;
    }

    BellGenerator generator(material, BellMeshConfig{15});
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) {
        std::cerr << "failed baseline mesh\n";
        return 2;
    }
    const auto& mesh = generator.getMeshData();

    const std::vector<double> q_values{1.00, 0.95, 0.90, 0.85, 0.80, 0.70};
    std::vector<std::vector<double>> frequencies(4);
    std::cout << std::setprecision(12);
    std::cout << "ellipticity_sweep_definition: q=b/a, y_physical=q*y_computational\n";
    std::cout << "mode,q,frequency_hz,relative_to_q1,mass_ratio,energy_ratio\n";

    for (size_t m = 0; m < 4; ++m) {
        frequencies[m].reserve(q_values.size());
        for (double q : q_values) {
            const BellOptimizer::Math::CoordinateScale3D inverse_map{1.0, q, 1.0};
            const auto r = BellOptimizer::Math::evaluate_mechanical_morphing(
                mesh, material, base.modes[m], inverse_map);
            assert(std::isfinite(r.morphed_frequency_hz) && r.morphed_frequency_hz > 0.0);
            assert(std::isfinite(r.morphed_mass) && r.morphed_mass > 0.0);
            assert(std::isfinite(r.morphed_strain_energy) && r.morphed_strain_energy > 0.0);
            frequencies[m].push_back(r.morphed_frequency_hz);
            const double f0 = frequencies[m].front();
            std::cout << (m + 1) << ',' << q << ',' << r.morphed_frequency_hz << ','
                      << (r.morphed_frequency_hz / f0 - 1.0) << ','
                      << (r.morphed_mass / r.base_mass) << ','
                      << (r.morphed_strain_energy / r.base_strain_energy) << '\n';
        }
    }

    // q=1 is the identity gate. The remaining sweep is deliberately diagnostic:
    // without an independent non-axisymmetric reference it must not be presented
    // as an accuracy claim.
    for (size_t m = 0; m < 4; ++m) {
        assert(std::abs(frequencies[m][0] - base.modes[m].frequency_hz) /
               std::max(base.modes[m].frequency_hz, 1e-12) < 2e-7);
        for (double f : frequencies[m]) assert(std::isfinite(f) && f > 0.0);
    }

    std::cout << "ellipticity_sweep_contract: PASS\n";
    return 0;
}
