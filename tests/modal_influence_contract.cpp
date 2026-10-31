#include "bell_generator.h"
#include "fem/internal_modal_solver.hpp"
#include "excitation/meridional_influence.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
class UniformRadialInfluence final : public BellOptimizer::ModalInfluenceProvider {
public:
    BellOptimizer::ModalInfluence make(const BellMeshData& mesh) const override {
        BellOptimizer::ModalInfluence x;
        x.id = name();
        x.radial.assign(mesh.nodes.size(), 1.0);
        x.axial.assign(mesh.nodes.size(), 0.0);
        return x;
    }
    const char* name() const noexcept override { return "uniform_radial_test"; }
};

BellOptimizer::ModalData solve(const BellGeometry& geometry,
                               const BellMaterial& material,
                               const BellOptimizer::ModalInfluenceProvider& provider) {
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    assert(generator.buildMesh());
    BellOptimizer::InternalModalSolver solver(true, 6);
    BellOptimizer::ModalData result;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "modal_influence", &provider};
    assert(solver.solve(problem, result));
    return result;
}
}

int main() {
    BellMaterial material;
    BellGeometry geometry{1000.0, 500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0}};

    // There is no hidden/default compatibility provider anymore. The current
    // production mathematical influence is explicitly selected by the caller.
    const BellOptimizer::Excitation::MeridionalInfluence meridional;
    const auto meridional_data = solve(geometry, material, meridional);

    // A different neutral influence changes only generalized excitation, not
    // the FEM eigenfrequencies. Uniform radial influence equals the part_x
    // projection by construction.
    const UniformRadialInfluence radial;
    const auto radial_data = solve(geometry, material, radial);
    assert(radial_data.modes.size() == meridional_data.modes.size());
    for (size_t i = 0; i < radial_data.modes.size(); ++i) {
        assert(std::abs(radial_data.modes[i].frequency_hz - meridional_data.modes[i].frequency_hz) < 1e-10);
        assert(std::abs(radial_data.modes[i].part_rz - radial_data.modes[i].part_x) < 1e-12);
        assert(std::abs(radial_data.modes[i].mass_rz - radial_data.modes[i].part_x * radial_data.modes[i].part_x) < 1e-12);
    }

    std::cout << "modal influence contract: PASS\n";
    return 0;
}
