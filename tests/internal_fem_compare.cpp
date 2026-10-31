#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "ccx/ccx_parser.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#ifndef BELL_SOURCE_DIR
#define BELL_SOURCE_DIR "."
#endif

int main(int argc, char** argv) {
    const std::string dat = argc > 1 ? argv[1] : std::string(BELL_SOURCE_DIR) + "/tests/reference/ccx/best_bell.dat";

    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return 2;

    BellOptimizer::InternalModalSolver internal(false, 10);
    BellOptimizer::ModalData a;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "validation", &modal_influence_1};
    if (!internal.solve(problem, a)) return 2;

    if (argc == 1) {
        for (const auto& mode : a.modes) {
            if (!(mode.frequency_hz > 0.0) || !std::isfinite(mode.frequency_hz) ||
                !std::isfinite(mode.mass_rz) || !std::isfinite(mode.acoustic_volume_displacement))
                return 5;
        }
        std::cout << "internal FEM sanity: " << a.modes.size()
                  << " finite positive modes on corrected normal-offset mesh\n";
        return a.modes.size() >= 5 ? 0 : 6;
    }

    const auto b = BellOptimizer::CCXParser::parse_file(dat, true);
    const size_t n = std::min(a.modes.size(), b.modes.size());
    std::cout << std::setprecision(10);
    std::cout << "mode,internal_hz,ccx_hz,relative_error\n";
    double maxerr = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double e = std::abs(a.modes[i].frequency_hz - b.modes[i].frequency_hz) / b.modes[i].frequency_hz;
        maxerr = std::max(maxerr, e);
        std::cout << (i + 1) << ',' << a.modes[i].frequency_hz << ',' << b.modes[i].frequency_hz << ',' << e << '\n';
    }
    std::cout << "max_relative_error=" << maxerr << "\n";
    if (n >= 1) {
        std::cout << "mode,internal_mass_rz,ccx_mass_rz,internal_amplitude,ccx_amplitude,amplitude_error\n";
        double max_amp_error = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const double e = std::abs(a.modes[i].amplitude - b.modes[i].amplitude);
            max_amp_error = std::max(max_amp_error, e);
            std::cout << (i + 1) << ',' << a.modes[i].mass_rz << ',' << b.modes[i].mass_rz << ','
                      << a.modes[i].amplitude << ',' << b.modes[i].amplitude << ',' << e << '\n';
        }
        std::cout << "max_amplitude_abs_error=" << max_amp_error << "\n";
        if (max_amp_error > 2.0e-3) return 4;
    }
    return (n >= 5 && maxerr < 1.0e-3) ? 0 : 3;
}
