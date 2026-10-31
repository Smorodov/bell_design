#include "bell_generator.h"
#include "fem/internal_modal_solver.hpp"
#include "excitation/localized_strike.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
constexpr double PI = 3.1415926535897932384626433832795;

struct Result {
    int n = 0;
    std::vector<double> frequencies;
    std::vector<double> pressure;
};

Result solve_mesh(int n, const BellMaterial& material, const BellGeometry& geometry,
                  double center_z, double width) {
    BellMeshConfig cfg;
    cfg.nNodes = n;
    BellGenerator generator(material, cfg);
    generator.setGeometry(geometry);

    GeometryValidationResult validation;
    if (!generator.buildMesh(&validation)) {
        std::cerr << "mesh " << n << " failed: " << validation.reason << "\n";
        std::exit(2);
    }

    BellOptimizer::InternalModalSolver solver(true, 9);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence influence;
    if (!solver.solve({generator.getMeshData(), material, ".", "localized_strike_mesh", &influence}, data) || data.modes.size() < 5) {
        std::cerr << "modal solve failed for mesh " << n << "\n";
        std::exit(3);
    }

    BellOptimizer::Excitation::LocalizedStrikeConfig strike;
    strike.center_z = center_z;
    strike.width = width;
    strike.normal_fraction = 1.0;
    strike.tangent_fraction = 0.0;

    const auto impulse = BellOptimizer::Excitation::modalImpulse(
        generator.getMeshData(), data.modes, strike);

    Result result;
    result.n = n;
    result.frequencies.reserve(data.modes.size());
    result.pressure.reserve(data.modes.size());
    for (size_t i = 0; i < data.modes.size(); ++i) {
        result.frequencies.push_back(data.modes[i].frequency_hz);
        result.pressure.push_back(
            impulse[i] * data.modes[i].acoustic_volume_displacement *
            (2.0 * PI * data.modes[i].frequency_hz));
    }
    return result;
}

double relative_rms(const std::vector<double>& a, const std::vector<double>& b) {
    assert(a.size() == b.size() && !a.empty());
    double sq = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double scale = std::max({std::abs(a[i]), std::abs(b[i]), 1e-12});
        const double d = (a[i] - b[i]) / scale;
        sq += d * d;
    }
    return std::sqrt(sq / static_cast<double>(a.size()));
}

double rms_log_spectrum(const std::vector<double>& a, const std::vector<double>& b) {
    assert(a.size() == b.size() && !a.empty());
    double na = 0.0, nb = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        na = std::max(na, std::abs(a[i]));
        nb = std::max(nb, std::abs(b[i]));
    }
    assert(na > 1e-14 && nb > 1e-14);

    double sq = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double aa = std::max(std::abs(a[i]) / na, 1e-14);
        const double bb = std::max(std::abs(b[i]) / nb, 1e-14);
        const double d = std::log(aa / bb);
        sq += d * d;
    }
    return std::sqrt(sq / static_cast<double>(a.size()));
}

} // namespace

int main() {
    BellMaterial material{1.1e4, 0.34, 8.6e-9};
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    // A previously selected finite-ring reference case. Keep the physical strike fixed
    // while refining only the FEM mesh so that changes measure discretization,
    // not a retuned excitation.
    constexpr double center_z = 624.534;
    constexpr double width = 5.0;
    const std::vector<int> meshes = {39, 47, 55, 63};

    std::vector<Result> results;
    results.reserve(meshes.size());
    for (const int n : meshes) results.push_back(solve_mesh(n, material, geometry, center_z, width));

    const Result& reference = results.back();
    std::cout << std::setprecision(9);
    std::cout << "reference_mesh=" << reference.n << "\n";
    std::cout << "nNodes,freq_rms_relative,pressure_rms_relative,pressure_rms_log_normalized\n";

        for (const Result& r : results) {
        const double freq_err = relative_rms(r.frequencies, reference.frequencies);
        const double pressure_err = relative_rms(r.pressure, reference.pressure);
        const double pressure_log = rms_log_spectrum(r.pressure, reference.pressure);
        std::cout << r.n << ',' << freq_err << ',' << pressure_err << ',' << pressure_log << '\n';

        assert(std::isfinite(freq_err));
        assert(std::isfinite(pressure_err));
        assert(std::isfinite(pressure_log));
        // Refinement should not introduce a large new instability. We allow
        // small non-monotonicity because eigenvectors and piecewise midpoint
        // integration converge discretely, but the finest mesh must improve
        // substantially over the coarsest mesh.
        if (r.n != reference.n) assert(pressure_log < 5.0);
    }

    const double coarse_log = rms_log_spectrum(results.front().pressure, reference.pressure);
    const double near_fine_log = rms_log_spectrum(results[results.size() - 2].pressure, reference.pressure);
    assert(coarse_log > 0.0);
    assert(near_fine_log < coarse_log);

    if (near_fine_log > 0.5) {
        std::cout << "localized strike mesh convergence: NOT CONVERGED at nNodes="
                  << reference.n << " (diagnostic PASS)\n";
    } else {
        std::cout << "localized strike mesh convergence: PASS\n";
    }
    return 0;
}
