#include "core/bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace BellOptimizer;

static double reference_profile(const char* family, double x) {
    // Smooth dimensionless reference profiles. They are deliberately denser
    // than the optimizer parameterization, so this test answers whether the
    // chosen knot count can represent a smooth forming profile rather than
    // merely reproducing another five-knot spline.
    if (std::string(family) == "european") {
        return x + 0.12 * x * (1.0 - x) * (1.0 - 2.0 * x);
    }
    if (std::string(family) == "chinese") {
        return x + 0.16 * x * (1.0 - x) * (1.0 - 1.5 * x);
    }
    return x + 0.10 * x * (1.0 - x) * (1.0 - 0.5 * x);
}

static BellGeometry make_base(const char* family, int knot_count) {
    BellGeometry g;
    g.H = 1000.0;
    g.R_mouth = (std::string(family) == "chinese") ? 560.0 : 500.0;
    g.t_top = (std::string(family) == "chinese") ? 11.0 : 10.0;
    g.t_bottom = (std::string(family) == "chinese") ? 32.0 : 30.0;
    g.r_knots.reserve(static_cast<size_t>(knot_count));
    for (int i = 0; i < knot_count; ++i) {
        const double x = static_cast<double>(i) / static_cast<double>(knot_count - 1);
        g.r_knots.push_back(reference_profile(family, x));
    }
    return g;
}

static std::vector<double> solve_freqs(const BellMaterial& material, const BellGeometry& g) {
    BellMeshConfig cfg;
    cfg.nNodes = 41;
    BellGenerator gen(material, cfg);
    gen.setGeometry(g);
    GeometryValidationResult validation;
    if (!gen.buildMesh(&validation)) {
        std::cerr << "mesh failure: " << validation.reason << "\n";
        return {};
    }
    InternalModalSolver solver;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    ModalProblem problem{gen.getMeshData(), material, ".", "spline_convergence", &modal_influence_1};
    ModalData data;
    if (!solver.solve(problem, data)) return {};
    std::vector<double> f;
    for (size_t i = 0; i < std::min<size_t>(5, data.modes.size()); ++i)
        f.push_back(data.modes[i].frequency_hz);
    return f;
}

int main() {
    BellMaterial material;
    BellMeshConfig cfg;
    cfg.nNodes = 51;

    double worst_profile_deviation = 0.0;
    double worst_frequency_error = 0.0;
    std::cout << std::setprecision(12);
    std::cout << "family, knots, max_profile_relative_error, max_frequency_relative_error\n";

    for (const char* family : {"european", "chinese", "neutral"}) {
        // 41 knots are used as a dense spline reference to keep the comparison
        // finite-dimensional while resolving the smooth analytic prior closely.
        const BellGeometry reference_geometry = make_base(family, 41);
        const BellGenerator reference_generator(material, cfg);
        const_cast<BellGenerator&>(reference_generator).setGeometry(reference_geometry);
        GeometryValidationResult validation;
        if (!const_cast<BellGenerator&>(reference_generator).buildMesh(&validation)) {
            std::cerr << family << " reference mesh failure: " << validation.reason << "\n";
            return 1;
        }
        const auto reference_freqs = solve_freqs(material, reference_geometry);
        if (reference_freqs.size() < 5) return 2;

        for (int knots : {5, 7, 9, 11, 13, 17}) {
            const BellGeometry candidate = make_base(family, knots);
            BellGenerator check(material, cfg);
            check.setGeometry(candidate);
            if (!check.buildMesh(&validation)) {
                std::cerr << family << " candidate mesh failure at " << knots
                          << " knots: " << validation.reason << "\n";
                return 3;
            }

            double profile_error = 0.0;
            for (int i = 1; i < 10000; ++i) {
                const double z = reference_geometry.H * static_cast<double>(i) / 10000.0;
                const double x = z / reference_geometry.H;
                const double ref = reference_profile(family, x) * reference_geometry.R_mouth;
                const double got = check.getOuterRadius(z);
                profile_error = std::max(profile_error,
                    std::abs(got - ref) / std::max(1.0, std::abs(ref)));
            }

            const auto f = solve_freqs(material, candidate);
            if (f.size() < 5) return 4;
            double freq_error = 0.0;
            for (size_t i = 0; i < 5; ++i) {
                freq_error = std::max(freq_error,
                    std::abs(f[i] - reference_freqs[i]) /
                    std::max(1.0, std::abs(reference_freqs[i])));
            }

            worst_profile_deviation = std::max(worst_profile_deviation, profile_error);
            worst_frequency_error = std::max(worst_frequency_error, freq_error);
            std::cout << family << ", " << knots << ", " << profile_error
                      << ", " << freq_error << "\n";
        }
    }

    // The audit establishes that five knots are numerically adequate for the
    // tested low modes, while the optimizer uses seven knots as a conservative
    // production parameterization to retain additional shape freedom.
    if (!(worst_frequency_error < 5e-3)) {
        std::cerr << "Spline convergence audit failed: increase profile knots.\n";
        return 5;
    }
    std::cout << "PASS: current smooth profile families are represented within the "
              << "0.5% first-five-mode frequency criterion.\n";
    std::cout << "Worst profile relative deviation = " << worst_profile_deviation << "\n";
    std::cout << "Worst first-five-mode frequency deviation = " << worst_frequency_error << "\n";
    return 0;
}
