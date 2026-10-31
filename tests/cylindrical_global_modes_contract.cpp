#include "fem/cylindrical_global_modes.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

static BellMeshData make_patch(int nr, int nz) {
    BellMeshData mesh;
    for (int j = 0; j <= nz; ++j) {
        const double z = 0.10 * static_cast<double>(j) / static_cast<double>(nz);
        for (int i = 0; i <= nr; ++i) {
            const double r = 0.20 + 0.10 * static_cast<double>(i) / static_cast<double>(nr);
            mesh.nodes.push_back({static_cast<int>(mesh.nodes.size()) + 1, r, z});
        }
    }
    const int stride = nr + 1;
    for (int j = 0; j < nz; ++j) for (int i = 0; i < nr; ++i) {
        const int n1 = j * stride + i + 1;
        const int n2 = n1 + 1;
        const int n4 = (j + 1) * stride + i + 1;
        const int n3 = n4 + 1;
        mesh.elements.push_back({n1, n2, n3, n4});
    }
    return mesh;
}

static std::vector<int> clamp_lower_edge(const BellMeshData& mesh, std::size_t dofs_per_node) {
    std::vector<int> fixed;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (std::abs(mesh.nodes[i].z) < 1e-12)
            for (std::size_t d = 0; d < dofs_per_node; ++d)
                fixed.push_back(static_cast<int>(i * dofs_per_node + d));
    }
    return fixed;
}

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::FEM;
    constexpr double E = 1.0e11, nu = 0.30, rho = 7800.0;
    const auto coarse_mesh = make_patch(2, 2);
    const auto fine_mesh = make_patch(4, 4);
    const auto coarse = assemble_cylindrical_global_system(coarse_mesh, E, nu, rho, {0}, 32);
    const auto fine = assemble_cylindrical_global_system(fine_mesh, E, nu, rho, {0}, 32);
    const auto circle_explicit = assemble_cylindrical_global_system(coarse_mesh, E, nu, rho, {0}, 32, 0.0, 1.0);
    const auto circle_multi = assemble_cylindrical_global_system(coarse_mesh, E, nu, rho, {0,1,2}, 32, 0.0, 1.0);
    const auto ellipse32 = assemble_cylindrical_global_system(coarse_mesh, E, nu, rho, {0,1,2}, 32, 0.0, 0.85);
    const auto ellipse64 = assemble_cylindrical_global_system(coarse_mesh, E, nu, rho, {0,1,2}, 64, 0.0, 0.85);
    auto relative_matrix_difference = [](const detail::SparseMatrix& a, const detail::SparseMatrix& b) {
        double diff2 = 0.0, scale2 = 0.0;
        for (int i = 0; i < a.n; ++i) for (int j = 0; j < a.n; ++j) {
            const double av = a.get(i,j), bv = b.get(i,j);
            diff2 += (av-bv)*(av-bv); scale2 += bv*bv;
        }
        return std::sqrt(diff2) / std::max(std::sqrt(scale2), 1e-300);
    };
    require(relative_matrix_difference(coarse.stiffness, circle_explicit.stiffness) == 0.0 &&
            relative_matrix_difference(coarse.mass, circle_explicit.mass) == 0.0,
            "explicit circular limit must reproduce the default assembly exactly");
    const double ellipse_k_quadrature_error = relative_matrix_difference(ellipse32.stiffness, ellipse64.stiffness);
    const double ellipse_m_quadrature_error = relative_matrix_difference(ellipse32.mass, ellipse64.mass);
    require(ellipse_k_quadrature_error < 1e-3 && ellipse_m_quadrature_error < 1e-3,
            "elliptical angular integration did not converge between 32 and 64 samples");
    require(relative_matrix_difference(circle_multi.stiffness, ellipse32.stiffness) > 1e-3,
            "ellipticity should change the assembled stiffness matrix");
    const auto circle_result = solve_cylindrical_global_modes(
        circle_multi, clamp_lower_edge(coarse_mesh, circle_multi.dofs_per_node), 2, 700, 5e-5, true);
    const auto ellipse_result = solve_cylindrical_global_modes(
        ellipse32, clamp_lower_edge(coarse_mesh, ellipse32.dofs_per_node), 2, 700, 5e-5, true);
    require(circle_result.converged && ellipse_result.converged &&
            circle_result.frequencies_hz.size() == 2 && ellipse_result.frequencies_hz.size() == 2,
            "circular/elliptical coupled-harmonic modal solve did not converge");
    const double first_frequency_change = std::abs(ellipse_result.frequencies_hz[0] - circle_result.frequencies_hz[0]) /
                                          circle_result.frequencies_hz[0];
    require(first_frequency_change > 1e-3,
            "ellipticity should change at least one low coupled-harmonic frequency measurably");
    const auto coarse_result = solve_cylindrical_global_modes(
        coarse, clamp_lower_edge(coarse_mesh, coarse.dofs_per_node), 2, 500, 5e-5, false);
    const auto fine_result = solve_cylindrical_global_modes(
        fine, clamp_lower_edge(fine_mesh, fine.dofs_per_node), 2, 700, 5e-5, true);
    require(coarse_result.converged && fine_result.converged,
            "multi-element constrained generalized eigensolver did not converge");
    require(coarse_result.frequencies_hz.size() == 2 && fine_result.frequencies_hz.size() == 2,
            "wrong number of modal results");
    require(coarse_result.relative_residual < 1e-4 && fine_result.relative_residual < 1e-4,
            "generalized eigenpair residual too large");
    for (const auto* result : {&coarse_result, &fine_result}) {
        for (std::size_t i = 0; i < result->frequencies_hz.size(); ++i) {
            require(std::isfinite(result->frequencies_hz[i]) && result->frequencies_hz[i] > 0.0,
                    "frequency must be positive and finite");
            if (i) require(result->frequencies_hz[i] >= result->frequencies_hz[i-1],
                           "eigenfrequencies must be sorted");
            for (double v : result->mode_shapes[i]) require(std::isfinite(v), "mode shape must be finite");
        }
    }
    for (std::size_t mode = 0; mode < 2; ++mode) {
        for (int dof : clamp_lower_edge(coarse_mesh, coarse.dofs_per_node))
            require(coarse_result.mode_shapes[mode][static_cast<std::size_t>(dof)] == 0.0,
                    "coarse fixed DOF must remain exactly zero");
        for (int dof : clamp_lower_edge(fine_mesh, fine.dofs_per_node))
            require(fine_result.mode_shapes[mode][static_cast<std::size_t>(dof)] == 0.0,
                    "fine fixed DOF must remain exactly zero");
    }
    // Refinement sanity: the two lowest eigenfrequencies should not jump wildly
    // when the meridional mesh is refined. This is a coarse convergence gate, not
    // a claim that the 4x4 result is a fully converged bell frequency.
    for (std::size_t mode = 0; mode < 2; ++mode) {
        const double rel = std::abs(fine_result.frequencies_hz[mode] - coarse_result.frequencies_hz[mode]) /
                           fine_result.frequencies_hz[mode];
        require(rel < 0.35, "mesh refinement changed a low mode excessively");
    }
    // Circular symmetry produces degenerate sine/cosine mode pairs. Ellipticity
    // breaks rotational symmetry, so those pairs split; do not assume that every
    // individual sorted frequency varies monotonically with q.
    const auto circle_spectrum = solve_cylindrical_global_modes(
        circle_multi, clamp_lower_edge(coarse_mesh, circle_multi.dofs_per_node),
        6, 1200, 1e-4, true);
    require(circle_spectrum.converged && circle_spectrum.frequencies_hz.size() == 6,
            "circular low-spectrum solve failed");
    const double circular_pair_split = std::abs(circle_spectrum.frequencies_hz[1] -
                                                 circle_spectrum.frequencies_hz[0]) /
                                       circle_spectrum.frequencies_hz[0];
    require(circular_pair_split < 1e-3,
            "circular rotational symmetry did not preserve the lowest modal pair");
    for (double q : {0.95, 0.85, 0.70}) {
        const auto swept_system = assemble_cylindrical_global_system(
            coarse_mesh, E, nu, rho, {0,1,2}, 64, 0.0, q);
        const auto swept = solve_cylindrical_global_modes(
            swept_system, clamp_lower_edge(coarse_mesh, swept_system.dofs_per_node),
            6, 1200, 1e-4, true);
        require(swept.converged && swept.frequencies_hz.size() == 6,
                "ellipticity sweep modal solve failed");
        require(swept.relative_residual < 1e-4,
                "ellipticity sweep residual exceeded tolerance");
        for (std::size_t i = 0; i < swept.frequencies_hz.size(); ++i) {
            require(std::isfinite(swept.frequencies_hz[i]) && swept.frequencies_hz[i] > 0.0,
                    "ellipticity sweep returned invalid frequency");
            if (i) require(swept.frequencies_hz[i] >= swept.frequencies_hz[i-1],
                           "ellipticity sweep frequencies are not sorted");
        }
        const double pair_split = std::abs(swept.frequencies_hz[1] - swept.frequencies_hz[0]) /
                                  swept.frequencies_hz[0];
        require(pair_split > circular_pair_split + 5e-4,
                "ellipticity did not measurably split the circular lowest modal pair");
        std::cout << "ellipticity_q=" << q << " first_pair_Hz="
                  << swept.frequencies_hz[0] << ',' << swept.frequencies_hz[1]
                  << " relative_split=" << pair_split << " residual=" << swept.relative_residual << '\n';
    }
    bool duplicate_rejected = false;
    try { (void)solve_cylindrical_global_modes(coarse, {0,0}, 1); }
    catch (const std::invalid_argument&) { duplicate_rejected = true; }
    require(duplicate_rejected, "duplicate fixed DOF was not rejected");
    std::cout << "cylindrical_global_modes_contract: 2x2/4x4 refinement, constrained sparse modes PASS\n"
              << "coarse_Hz=" << coarse_result.frequencies_hz[0] << ',' << coarse_result.frequencies_hz[1]
              << " fine_Hz=" << fine_result.frequencies_hz[0] << ',' << fine_result.frequencies_hz[1]
              << " residuals=" << coarse_result.relative_residual << ',' << fine_result.relative_residual << '\n'
              << "circle_q100_Hz=" << circle_result.frequencies_hz[0] << ',' << circle_result.frequencies_hz[1] << '\n'
              << "ellipse_q085_Hz=" << ellipse_result.frequencies_hz[0] << ',' << ellipse_result.frequencies_hz[1]
              << " quadrature_errors_KM=" << ellipse_k_quadrature_error << ',' << ellipse_m_quadrature_error << '\n';
}
