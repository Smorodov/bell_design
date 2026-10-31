#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Case { const char* name; BellGeometry geometry; };

std::vector<Case> cases() {
    BellGeometry a{1000.0, 500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0}};
    BellGeometry b{1000.0, 600.0, 12.0, 35.0, {0.0, 0.12, 0.35, 0.65, 1.0}};
    BellGeometry c{1000.0, 400.0, 8.0, 25.0, {0.0, 0.18, 0.45, 0.75, 1.0}};
    BellGeometry d{1500.0, 750.0, 15.0, 45.0, {0.0, 0.14, 0.38, 0.68, 1.0}};
    BellGeometry e{600.0, 300.0, 6.0, 18.0, {0.0, 0.16, 0.42, 0.72, 1.0}};
    return {{"g1", a}, {"g2", b}, {"g3", c}, {"g4", d}, {"g5", e}};
}

bool solve(const BellGeometry& geometry, const BellMaterial& material, int nNodes,
           BellOptimizer::ModalData& result) {
    BellMeshConfig config;
    config.nNodes = nNodes;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return false;
    BellOptimizer::InternalModalSolver solver(true, 8);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "validation", &modal_influence_1};
    return solver.solve(problem, result);
}

double rel(double a, double b) {
    return std::abs(a - b) / std::max(std::abs(b), 1e-12);
}

} // namespace

int main() {
    BellMaterial base_material;
    const auto geometries = cases();
    const double eps_scale = 2.0e-3;
    const double eps_invariance = 2.0e-3;
    const double eps_refinement = 1.5e-1;
    double worst_scale = 0.0;
    double worst_ratio = 0.0;
    double worst_refinement = 0.0;

    std::cout << std::setprecision(10);
    std::cout << "case,scale,mode,frequency_hz,expected_hz,relative_error\n";

    for (const auto& item : geometries) {
        BellOptimizer::ModalData base;
        if (!solve(item.geometry, base_material, 15, base) || base.modes.size() < 5) {
            std::cerr << "failed baseline: " << item.name << '\n';
            return 2;
        }

        const double base_prime = base.modes[1].frequency_hz;
        for (double scale : {0.75, 1.25}) {
            BellGeometry scaled = item.geometry;
            scaled.scale(scale);
            BellOptimizer::ModalData data;
            if (!solve(scaled, base_material, 15, data) || data.modes.size() < 5) {
                std::cerr << "failed scaled case: " << item.name << '\n';
                return 2;
            }
            for (size_t i = 0; i < 5; ++i) {
                const double expected = base.modes[i].frequency_hz / scale;
                const double error = rel(data.modes[i].frequency_hz, expected);
                worst_scale = std::max(worst_scale, error);
                std::cout << item.name << ',' << scale << ',' << (i + 1) << ','
                          << data.modes[i].frequency_hz << ',' << expected << ',' << error << '\n';
            }
            for (size_t i = 0; i < 5; ++i) {
                const double r0 = base.modes[i].frequency_hz / base_prime;
                const double r1 = data.modes[i].frequency_hz / data.modes[1].frequency_hz;
                worst_ratio = std::max(worst_ratio, std::abs(r1 - r0));
            }
        }

        if (item.name == std::string("g1")) {
            for (double factor : {0.5, 2.0}) {
                BellMaterial m = base_material;
                m.E *= factor;
                BellOptimizer::ModalData data;
                if (!solve(item.geometry, m, 15, data) || data.modes.size() < 5) return 2;
                for (size_t i = 0; i < 5; ++i) {
                    const double expected = base.modes[i].frequency_hz * std::sqrt(factor);
                    worst_scale = std::max(worst_scale, rel(data.modes[i].frequency_hz, expected));
                }
                m = base_material;
                m.rho *= factor;
                if (!solve(item.geometry, m, 15, data) || data.modes.size() < 5) return 2;
                for (size_t i = 0; i < 5; ++i) {
                    const double expected = base.modes[i].frequency_hz / std::sqrt(factor);
                    worst_scale = std::max(worst_scale, rel(data.modes[i].frequency_hz, expected));
                }
            }
        }

        if (item.name == std::string("g1")) {
            BellOptimizer::ModalData coarse, fine, finer;
            if (!solve(item.geometry, base_material, 10, coarse) ||
                !solve(item.geometry, base_material, 15, fine) ||
                !solve(item.geometry, base_material, 20, finer) ||
                coarse.modes.size() < 5 || fine.modes.size() < 5 || finer.modes.size() < 5) return 2;
            for (size_t i = 0; i < 5; ++i) {
                const double coarse_to_fine = rel(fine.modes[i].frequency_hz, coarse.modes[i].frequency_hz);
                const double fine_to_finer = rel(finer.modes[i].frequency_hz, fine.modes[i].frequency_hz);
                worst_refinement = std::max(worst_refinement, fine_to_finer);
                std::cout << "refinement,g1," << (i + 1) << ','
                          << coarse_to_fine << ',' << fine_to_finer << '\n';
            }
        }

    }

    std::cout << "worst_scaling_relative_error=" << worst_scale << '\n';
    std::cout << "worst_ratio_invariance_error=" << worst_ratio << '\n';
    std::cout << "worst_refinement_relative_error=" << worst_refinement << '\n';

    if (worst_scale > eps_scale) return 3;
    if (worst_ratio > eps_invariance) return 4;
    if (worst_refinement > eps_refinement) return 5;
    return 0;
}
