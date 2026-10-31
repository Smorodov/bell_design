#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

static bool valid_result(const BellOptimizer::ModalData& result) {
    if (!result.matrix_audit.available || result.modes.size() < 5) return false;
    for (const auto& mode : result.modes) {
        if (!(mode.frequency_hz > 0.0) || !std::isfinite(mode.frequency_hz)) return false;
        if (!std::isfinite(mode.eigen_residual_relative) || mode.eigen_residual_relative > 1e-5) return false;
        if (!std::isfinite(mode.modal_mass) || std::abs(mode.modal_mass - 1.0) > 1e-7) return false;
    }
    return true;
}

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 121;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};

    auto run = [&](bool rcm, BellOptimizer::ModalData& out) {
        BellGenerator generator(material, config);
        generator.setGeometry(geometry);
        if (!generator.buildMesh()) return false;
        BellOptimizer::InternalModalSolver solver(true, 5, true, 128, rcm);
        BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
        BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", rcm ? "rcm_on" : "rcm_off", &modal_influence_1};
        const auto begin = std::chrono::steady_clock::now();
        const bool ok = solver.solve(problem, out);
        const auto end = std::chrono::steady_clock::now();
        if (!ok) return false;
        std::cout << (rcm ? "rcm_on" : "rcm_off")
                  << "_time_ms="
                  << std::chrono::duration<double, std::milli>(end - begin).count() << '\n';
        return valid_result(out);
    };

    BellOptimizer::ModalData off, on;
    if (!run(false, off) || !run(true, on)) return 1;
    if (!on.matrix_audit.production_rcm_ordering) return 2;
    if (off.matrix_audit.production_rcm_ordering) return 3;
    if (off.matrix_audit.stiffness.ic0_factor_nonzeros == 0 || on.matrix_audit.stiffness.ic0_factor_nonzeros == 0) return 4;
    if (on.modes.size() != off.modes.size()) return 5;

    double max_frequency_error = 0.0;
    double min_mode_correlation = 1.0;
    for (std::size_t i = 0; i < off.modes.size(); ++i) {
        const double f0 = off.modes[i].frequency_hz;
        max_frequency_error = std::max(max_frequency_error, std::abs(on.modes[i].frequency_hz - f0) /
            std::max(std::abs(f0), 1e-12));
        if (off.modes[i].nodal_radial.size() == on.modes[i].nodal_radial.size() &&
            !off.modes[i].nodal_radial.empty()) {
            double dot = 0.0, a2 = 0.0, b2 = 0.0;
            for (std::size_t j = 0; j < off.modes[i].nodal_radial.size(); ++j) {
                const double a = off.modes[i].nodal_radial[j];
                const double b = on.modes[i].nodal_radial[j];
                dot += a*b; a2 += a*a; b2 += b*b;
            }
            if (a2 > 0.0 && b2 > 0.0)
                min_mode_correlation = std::min(min_mode_correlation, std::abs(dot) / std::sqrt(a2*b2));
        }
    }
    if (max_frequency_error > 1e-7 || min_mode_correlation < 0.999) return 6;

    std::cout << "rcm_solver_benchmark=PASS\n"
              << "reduced_dof=" << off.matrix_audit.reduced_dof << '\n'
              << "K_nnz=" << off.matrix_audit.stiffness.exact_nonzeros << '\n'
              << "M_nnz=" << off.matrix_audit.mass.exact_nonzeros << '\n'
              << "IC0_nnz_off=" << off.matrix_audit.stiffness.ic0_factor_nonzeros << '\n'
              << "IC0_nnz_on=" << on.matrix_audit.stiffness.ic0_factor_nonzeros << '\n'
              << "outer_iterations_off=" << off.matrix_audit.sparse_solver_iterations << '\n'
              << "outer_iterations_on=" << on.matrix_audit.sparse_solver_iterations << '\n'
              << "pcg_iterations_off=" << off.matrix_audit.sparse_pcg_iterations << '\n'
              << "pcg_iterations_on=" << on.matrix_audit.sparse_pcg_iterations << '\n'
              << "IC0_reduction_fraction="
              << (1.0 - static_cast<double>(on.matrix_audit.stiffness.ic0_factor_nonzeros) /
                  static_cast<double>(off.matrix_audit.stiffness.ic0_factor_nonzeros)) << '\n'
              << "max_frequency_relative_error=" << max_frequency_error << '\n'
              << "min_mode_correlation=" << min_mode_correlation << '\n';
    return 0;
}
