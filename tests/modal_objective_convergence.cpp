#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

struct Metrics {
    double prime = 0.0;
    std::array<double,4> freq_ratios{};
    std::array<double,4> amp_ratios{};
    double objective = 0.0;
};

bool solve(const BellGeometry& geometry, int nNodes, BellOptimizer::ModalData& result) {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = nNodes;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return false;
    BellOptimizer::InternalModalSolver solver(true, 8);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "objective_validation", &modal_influence_1};
    return solver.solve(problem, result);
}

double log_relative_sq(double value, double target) {
    if (!(value > 0.0) || !(target > 0.0) || !std::isfinite(value) || !std::isfinite(target)) {
        return 1e9;
    }
    const double e = std::log(value / target);
    return e * e;
}


Metrics measure(const BellOptimizer::ModalData& data) {
    Metrics m;
    if (data.modes.size() < 5) return m;

    const double prime_raw = data.modes[1].frequency_hz;
    const double scale = prime_raw / 440.0;
    m.prime = prime_raw / scale;

    const int freq_indices[4] = {0, 2, 3, 4};
    const double freq_targets[4] = {0.5, 1.25, 1.5, 2.0};
    const double amp_targets[4] = {0.8, 0.6, 0.4, 0.2};

    double freq_error = 0.0;
    for (int i = 0; i < 4; ++i) {
        m.freq_ratios[i] = data.modes[freq_indices[i]].frequency_hz / prime_raw;
        freq_error += log_relative_sq(m.freq_ratios[i], freq_targets[i]);
    }

    double amp_error = 0.0;
    for (int i = 0; i < 4; ++i) {
        m.amp_ratios[i] = data.modes[i + 1].amplitude / data.modes[0].amplitude;
        amp_error += log_relative_sq(m.amp_ratios[i], amp_targets[i]);
    }

    m.objective = freq_error + amp_error;
    return m;
}

double max_physical_metric_difference(const Metrics& a, const Metrics& b) {
    double d = std::abs(a.prime - b.prime);
    for (size_t i = 0; i < 4; ++i) {
        d = std::max(d, std::abs(a.freq_ratios[i] - b.freq_ratios[i]));
        d = std::max(d, std::abs(a.amp_ratios[i] - b.amp_ratios[i]));
    }
    return d;
}

} // namespace

int main() {
    BellGeometry geometry{1000.0, 500.0, 10.0, 30.0, {0.0, 0.15, 0.40, 0.70, 1.0}};
    const std::array<int,5> meshes = {10, 15, 20, 25, 30};
    std::array<Metrics,5> metrics{};

    std::cout << std::setprecision(10);
    std::cout << "nNodes,prime_hz,hum_prime,tertia_prime,quinta_prime,nominal_prime,";
    std::cout << "prime_hum,tertia_hum,quinta_hum,nominal_hum,objective\n";

    for (size_t k = 0; k < meshes.size(); ++k) {
        BellOptimizer::ModalData data;
        if (!solve(geometry, meshes[k], data)) {
            std::cerr << "solver failed for nNodes=" << meshes[k] << '\n';
            return 2;
        }
        metrics[k] = measure(data);
        const auto& m = metrics[k];
        std::cout << meshes[k] << ',' << m.prime;
        for (double x : m.freq_ratios) std::cout << ',' << x;
        for (double x : m.amp_ratios) std::cout << ',' << x;
        std::cout << ',' << m.objective << '\n';
    }

    const Metrics& reference = metrics.back();
    double worst = 0.0;
    for (size_t k = 0; k + 1 < metrics.size(); ++k)
        worst = std::max(worst, max_physical_metric_difference(metrics[k], reference));

    const double d10_15 = max_physical_metric_difference(metrics[0], metrics[1]);
    const double d15_20 = max_physical_metric_difference(metrics[1], metrics[2]);
    const double d20_25 = max_physical_metric_difference(metrics[2], metrics[3]);
    const double d25_30 = max_physical_metric_difference(metrics[3], metrics[4]);

    std::cout << "worst_difference_to_30=" << worst << '\n';
    std::cout << "difference_10_to_15=" << d10_15 << '\n';
    std::cout << "difference_15_to_20=" << d15_20 << '\n';
    std::cout << "difference_20_to_25=" << d20_25 << '\n';
    std::cout << "difference_25_to_30=" << d25_30 << '\n';

    // This is a convergence diagnostic, not a monotonicity requirement.
    // Correct normal-offset geometry and axis-focused node redistribution
    // change the sampling sequence; individual refinement steps may therefore
    // be non-monotone. The coarse-to-fine discrepancy must still shrink.
    if (!(d25_30 < d10_15 && d25_30 < 0.5 * d10_15)) return 3;
    return 0;
}
