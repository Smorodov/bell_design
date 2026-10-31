#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "fem/sparse_eigensolver.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

bool buildLayeredMesh(const BellMeshData& base, int profile_nodes, int layers, BellMeshData& out) {
    if (layers < 1 || base.outer_profile.size() != static_cast<std::size_t>(profile_nodes) ||
        base.inner_profile.size() != static_cast<std::size_t>(profile_nodes)) return false;
    out.outer_profile = base.outer_profile;
    out.inner_profile = base.inner_profile;
    out.total_arc_length = base.total_arc_length;
    out.min_thickness = base.min_thickness;
    out.max_thickness = base.max_thickness;
    out.nodes.reserve(static_cast<std::size_t>(profile_nodes) * static_cast<std::size_t>(layers + 1));
    out.elements.reserve(static_cast<std::size_t>(profile_nodes - 1) * static_cast<std::size_t>(layers));
    for (int layer = 0; layer <= layers; ++layer) {
        const double a = static_cast<double>(layer) / static_cast<double>(layers);
        for (int i = 0; i < profile_nodes; ++i) {
            const auto& o = base.outer_profile[static_cast<std::size_t>(i)];
            const auto& in = base.inner_profile[static_cast<std::size_t>(i)];
            out.nodes.push_back({static_cast<int>(out.nodes.size()) + 1,
                o.r + a * (in.r - o.r), o.z + a * (in.z - o.z)});
        }
    }
    for (int layer = 0; layer < layers; ++layer) {
        const int outer = layer * profile_nodes;
        const int inner = (layer + 1) * profile_nodes;
        for (int i = 0; i < profile_nodes - 1; ++i)
            out.elements.push_back({inner + i + 1, outer + i + 1, outer + i + 2, inner + i + 2});
    }
    return true;
}

bool clusteredSpectrumTest() {
    using BellOptimizer::detail::SparseMatrix;
    using BellOptimizer::detail::SparseTriplet;
    using BellOptimizer::detail::solve_lowest_sparse_modes;
    constexpr int n = 32;
    std::vector<SparseTriplet> kt, mt;
    for (int i = 0; i < n; ++i) {
        double lambda = 5.0 + static_cast<double>(i);
        if (i == 0) lambda = 1.0;
        if (i == 1) lambda = 1.0 + 1e-7;
        if (i == 2) lambda = 1.0 + 2e-7;
        kt.push_back({i, i, lambda});
        mt.push_back({i, i, 1.0});
    }
    const auto result = solve_lowest_sparse_modes(SparseMatrix::from_triplets(n, std::move(kt)),
        SparseMatrix::from_triplets(n, std::move(mt)), 3, 200, 1e-10);
    if (!result.converged || result.eigenvalues.size() != 3) { std::cerr << "cluster convergence=" << result.converged << " modes=" << result.eigenvalues.size() << " iterations=" << result.iterations << "\n"; return false; }
    const double expected[] = {1.0, 1.0 + 1e-7, 1.0 + 2e-7};
    for (int i = 0; i < 3; ++i) {
        if (std::abs(result.eigenvalues[static_cast<std::size_t>(i)] - expected[i]) > 2e-7) { std::cerr << "cluster eigenvalue " << i << "=" << result.eigenvalues[static_cast<std::size_t>(i)] << " expected=" << expected[i] << "\n"; return false; }
        std::vector<double> residual(static_cast<std::size_t>(n), 0.0);
        for (int j = 0; j < n; ++j) {
            double lambda = 5.0 + static_cast<double>(j);
            if (j == 0) lambda = 1.0;
            if (j == 1) lambda = 1.0 + 1e-7;
            if (j == 2) lambda = 1.0 + 2e-7;
            residual[static_cast<std::size_t>(j)] =
                (lambda - result.eigenvalues[static_cast<std::size_t>(i)]) *
                result.vectors[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
        double r2 = 0.0;
        for (double x : residual) r2 += x * x;
        if (std::sqrt(r2) > 2e-6) { std::cerr << "cluster residual mode=" << i << " residual=" << std::sqrt(r2) << " lambda=" << result.eigenvalues[static_cast<std::size_t>(i)] << "\n"; return false; }
    }
    if (!(result.max_mass_orthogonality_error < 1e-8) ||
        !(result.max_eigenpair_relative_residual < 1e-9) ||
        !std::isfinite(result.regularization_shift) || result.regularization_shift <= 0.0) {
        std::cerr << "cluster diagnostics invalid: orth=" << result.max_mass_orthogonality_error
                  << " residual=" << result.max_eigenpair_relative_residual
                  << " shift=" << result.regularization_shift << '\n';
        return false;
    }
    std::cout << "clustered_spectrum=PASS iterations=" << result.iterations
              << " orth=" << result.max_mass_orthogonality_error
              << " residual=" << result.max_eigenpair_relative_residual << '\n';
    return true;
}

bool runMeshCase(int profile_nodes, int layers, const BellMaterial& material) {
    const auto start = Clock::now();
    BellMeshConfig config;
    config.nNodes = profile_nodes;
    BellGenerator generator(material, config);
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
        22.85675347179014, 46.62524667849419,
        {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    BellMeshData mesh;
    if (!buildLayeredMesh(generator.getMeshData(), profile_nodes, layers, mesh)) return false;
    BellOptimizer::InternalModalSolver solver(true, 5, true, 0);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{mesh, material, ".", "sparse_solver_robustness", &modal_influence_1};
    BellOptimizer::ModalData result;
    if (!solver.solve(problem, result) || result.modes.size() != 5 || !result.matrix_audit.available) return false;
    double max_residual = 0.0;
    for (const auto& mode : result.modes) {
        if (!(mode.frequency_hz > 0.0) || !std::isfinite(mode.frequency_hz) ||
            !(mode.eigen_residual_relative < 2e-5) || std::abs(mode.modal_mass - 1.0) > 1e-7 ||
            !std::isfinite(mode.mass_orthogonality_error) || mode.mass_orthogonality_error > 2e-4) return false;
        max_residual = std::max(max_residual, mode.eigen_residual_relative);
    }
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start).count();
    const auto& ka = result.matrix_audit.stiffness;
    const auto& ma = result.matrix_audit.mass;
    std::cout << "mesh_case=PASS profile_nodes=" << profile_nodes << " layers=" << layers
              << " reduced_dof=" << result.matrix_audit.reduced_dof
              << " K_nnz=" << ka.exact_nonzeros << " M_nnz=" << ma.exact_nonzeros
              << " csr_bytes=" << (ka.estimated_csr_bytes + ma.estimated_csr_bytes)
              << " dense_matrix_pair_bytes=" << (ka.dense_bytes + ma.dense_bytes)
              << " f1_hz=" << result.modes.front().frequency_hz
              << " max_residual=" << max_residual << " elapsed_ms=" << ms << '\n';
    return true;
}
}

int main(int argc, char** argv) {
    BellMaterial material;
    if (argc == 3) {
        char* end_nodes = nullptr;
        char* end_layers = nullptr;
        const long nodes = std::strtol(argv[1], &end_nodes, 10);
        const long layers = std::strtol(argv[2], &end_layers, 10);
        if (!end_nodes || *end_nodes != '\0' || !end_layers || *end_layers != '\0' ||
            nodes < 3 || nodes > 100000 || layers < 1 || layers > 1000) {
            std::cerr << "usage: sparse_solver_robustness [profile_nodes layers]\n";
            return 64;
        }
        return runMeshCase(static_cast<int>(nodes), static_cast<int>(layers), material) ? 0 : 1;
    }
    if (argc != 1) {
        std::cerr << "usage: sparse_solver_robustness [profile_nodes layers]\n";
        return 64;
    }
    if (!clusteredSpectrumTest()) {
        std::cerr << "clustered-spectrum sparse eigensolver test failed\n";
        return 1;
    }
    for (int layers : {2, 3, 4})
        if (!runMeshCase(31, layers, material)) {
            std::cerr << "sparse robustness failed for profile_nodes=31 layers=" << layers << '\n';
            return 2;
        }
    if (!runMeshCase(61, 2, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=61 layers=2\n";
        return 3;
    }
    // Larger sparse-only cases exercise memory and convergence beyond the former test range.
    if (!runMeshCase(121, 2, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=121 layers=2\n";
        return 4;
    }
    if (!runMeshCase(121, 4, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=121 layers=4\n";
        return 5;
    }
    if (!runMeshCase(181, 2, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=181 layers=2\n";
        return 6;
    }
    if (!runMeshCase(181, 4, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=181 layers=4\n";
        return 7;
    }
    // 201 profile stations are the largest stable geometry-generated mesh in
    // this test geometry.  Larger station counts currently fail geometry
    // construction before FEM assembly, so they are not used as solver limits.
    if (!runMeshCase(201, 2, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=201 layers=2\n";
        return 8;
    }
    if (!runMeshCase(201, 4, material)) {
        std::cerr << "sparse robustness failed for profile_nodes=201 layers=4\n";
        return 9;
    }
    std::cout << "sparse_solver_robustness=PASS\n";
    return 0;
}
