#include "bell_generator.h"
#include "fem/internal_modal_solver.hpp"
#include "excitation/localized_strike.hpp"
#include "excitation/meridional_influence.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {
double norm(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}
double cosine(const std::vector<double>& a, const std::vector<double>& b) {
    double dot = 0.0, aa = 0.0, bb = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot += a[i] * b[i]; aa += a[i] * a[i]; bb += b[i] * b[i];
    }
    return dot / std::sqrt(aa * bb);
}
}

int main() {
    BellMaterial material{1.1e4, 0.34, 8.6e-9};
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellMeshConfig cfg; cfg.nNodes = 15;
    BellGenerator generator(material, cfg);
    generator.setGeometry(geometry);
    GeometryValidationResult validation;
    if (!generator.buildMesh(&validation)) return 2;

    BellOptimizer::InternalModalSolver solver(true, 7, true);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence influence;
    if (!solver.solve({generator.getMeshData(), material, ".", "localized_strike_direction", &influence}, data) || data.modes.size() < 5) return 3;

    BellOptimizer::Excitation::LocalizedStrikeConfig a;
    a.center_z = generator.getOuterProfile()[8].z;
    a.width = 20.0;
    a.normal_fraction = 1.0;
    a.tangent_fraction = 0.0;

    auto b = a;
    b.normal_fraction = 3.0;
    b.tangent_fraction = 4.0;

    auto c = a;
    c.normal_fraction = 0.6;
    c.tangent_fraction = 0.8;

    const auto ia = BellOptimizer::Excitation::modalImpulse(generator.getMeshData(), data.modes, a);
    const auto ib = BellOptimizer::Excitation::modalImpulse(generator.getMeshData(), data.modes, b);
    const auto ic = BellOptimizer::Excitation::modalImpulse(generator.getMeshData(), data.modes, c);

    if (norm(ia) < 1e-14 || norm(ib) < 1e-14 || norm(ic) < 1e-14) return 4;

    // (3,4) and (0.6,0.8) define exactly the same unit direction. Their modal
    // projections must therefore agree independently of the parameter scale.
    const double same_direction_error = 1.0 - cosine(ib, ic);
    if (!(same_direction_error < 1e-12)) {
        std::cerr << "same_direction_error=" << same_direction_error << "\n";
        return 5;
    }

    // A purely normal strike is a genuinely different excitation direction.
    const double different_direction_error = 1.0 - cosine(ia, ib);
    if (!(different_direction_error > 1e-6)) return 6;

    std::cout << "localized_strike_direction_contract: PASS\n";
    return 0;
}
