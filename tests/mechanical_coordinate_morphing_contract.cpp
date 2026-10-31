#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "math/modal_coordinate_morphing.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

double rel(double a, double b) {
    return std::abs(a - b) / std::max(std::abs(b), 1e-12);
}

bool solve(const BellGeometry& geometry, const BellMaterial& material,
           bool capture, BellOptimizer::ModalData& result) {
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return false;
    BellOptimizer::InternalModalSolver solver(true, 6, capture);
    BellOptimizer::Excitation::MeridionalInfluence influence;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "mechanical_morph", &influence};
    return solver.solve(problem, result);
}

} // namespace

int main() {
    BellMaterial material;
    BellGeometry geometry{1000.0, 500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0}};

    BellOptimizer::ModalData base;
    if (!solve(geometry, material, true, base) || base.modes.size() < 4) {
        std::cerr << "failed baseline solve, modes=" << base.modes.size() << "\n";
        return 2;
    }

    std::cout << std::setprecision(12);
    double worst_identity = 0.0;
    double worst_uniform = 0.0;

    // Identity must reproduce the FEM modal energy/frequency, not merely the
    // stored eigenvalue. This is the fundamental n=0 morphing gate.
    // Rebuild the mesh once so the morphing evaluator receives the exact mesh
    // used to generate the captured mode vectors.
    BellGenerator generator(material, BellMeshConfig{15});
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) { std::cerr << "failed mesh build\n"; return 2; }
    const auto& mesh = generator.getMeshData();

    for (size_t i = 0; i < 4; ++i) {
        const auto r = BellOptimizer::Math::evaluate_mechanical_morphing(
            mesh, material, base.modes[i], {});
        const double e_mass = rel(r.base_mass, 1.0);
        const double e_energy = rel(2.0 * r.base_strain_energy, base.modes[i].eigenvalue);
        const double e_freq = rel(r.base_frequency_hz, base.modes[i].frequency_hz);
        worst_identity = std::max(worst_identity, std::max({e_mass, e_energy, e_freq}));
        std::cout << "identity,mode=" << (i + 1)
                  << ",mass=" << r.base_mass
                  << ",energy=" << r.base_strain_energy
                  << ",frequency=" << r.base_frequency_hz
                  << ",reference=" << base.modes[i].frequency_hz << '\n';
    }

    // Uniform affine scaling has an exact frequency law f' = f/s. The reduced
    // modal calculation must reproduce the directly solved uniformly scaled
    // geometry to numerical accuracy.
    for (double s : {0.8, 1.2}) {
        BellOptimizer::Math::CoordinateScale3D scale{s, s, s};
        BellOptimizer::ModalData direct;
        BellGeometry scaled = geometry;
        scaled.scale(s);
        if (!solve(scaled, material, false, direct) || direct.modes.size() < 4) {
            std::cerr << "failed scaled solve\n"; return 2;
        }
        for (size_t i = 0; i < 4; ++i) {
            const auto r = BellOptimizer::Math::evaluate_mechanical_morphing(
                mesh, material, base.modes[i], scale);
            const double error = rel(r.morphed_frequency_hz, direct.modes[i].frequency_hz);
            worst_uniform = std::max(worst_uniform, error);
            const double expected = base.modes[i].frequency_hz / s;
            assert(rel(r.morphed_frequency_hz, expected) < 2e-7);
            std::cout << "uniform,scale=" << s << ",mode=" << (i + 1)
                      << ",morphed=" << r.morphed_frequency_hz
                      << ",direct=" << direct.modes[i].frequency_hz
                      << ",relative_error=" << error << '\n';
        }
    }

    // A transverse-only compression is deliberately not assigned a simple
    // frequency factor. It must produce a finite Rayleigh correction that is
    // mode-dependent; this is the intended reduced-order path toward elliptic
    // cross-sections.
    BellOptimizer::Math::CoordinateScale3D compressed_y{1.0, 0.85, 1.0};
    double first = 0.0;
    bool mode_dependence_seen = false;
    for (size_t i = 0; i < 4; ++i) {
        const auto r = BellOptimizer::Math::evaluate_mechanical_morphing(
            mesh, material, base.modes[i], compressed_y);
        assert(std::isfinite(r.morphed_frequency_hz) && r.morphed_frequency_hz > 0.0);
        if (i == 0) first = r.morphed_frequency_hz;
        else if (std::abs(r.morphed_frequency_hz - first) > 1e-6) mode_dependence_seen = true;
        std::cout << "transverse,scale_y=0.85,mode=" << (i + 1)
                  << ",frequency=" << r.morphed_frequency_hz
                  << ",relative_change=" << r.relative_frequency_change << '\n';
    }
    assert(mode_dependence_seen);

    std::cout << "worst_identity_relative_error=" << worst_identity << '\n';
    std::cout << "worst_uniform_relative_error=" << worst_uniform << '\n';
    assert(worst_identity < 2e-7);
    assert(worst_uniform < 2e-7);
    std::cout << "mechanical coordinate morphing contract: PASS\n";
    return 0;
}
