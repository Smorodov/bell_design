#include "fem/cylindrical_global_modes.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

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

static std::vector<int> clamp_lower_edge(const BellMeshData& mesh, std::size_t dpn) {
    std::vector<int> fixed;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (std::abs(mesh.nodes[i].z) < 1e-12)
            for (std::size_t d = 0; d < dpn; ++d)
                fixed.push_back(static_cast<int>(i * dpn + d));
    }
    return fixed;
}

static double max_relative_change(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double out = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        out = std::max(out, std::abs(a[i] - b[i]) / std::max(std::abs(b[i]), 1e-300));
    return out;
}


static double sparse_relative_difference(const detail::SparseMatrix& a,
                                         const detail::SparseMatrix& b) {
    require(a.n == b.n, "angular quadrature changed sparse matrix dimension");
    // Zero entries may be dropped at different quadrature resolutions, so the
    // CSR sparsity pattern itself is not a convergence contract. Merge the two
    // sorted row patterns and compare the logical matrix values in O(nnz).
    double diff2 = 0.0;
    double ref2 = 0.0;
    for (int row = 0; row < a.n; ++row) {
        std::size_t pa = a.row_offsets[static_cast<std::size_t>(row)];
        const std::size_t ea = a.row_offsets[static_cast<std::size_t>(row) + 1u];
        std::size_t pb = b.row_offsets[static_cast<std::size_t>(row)];
        const std::size_t eb = b.row_offsets[static_cast<std::size_t>(row) + 1u];
        while (pa < ea || pb < eb) {
            const int ca = pa < ea ? a.columns[pa] : std::numeric_limits<int>::max();
            const int cb = pb < eb ? b.columns[pb] : std::numeric_limits<int>::max();
            const int col = std::min(ca, cb);
            const double va = ca == col ? a.values[pa++] : 0.0;
            const double vb = cb == col ? b.values[pb++] : 0.0;
            const double d = va - vb;
            diff2 += d * d;
            ref2 += vb * vb;
        }
    }
    return std::sqrt(diff2 / std::max(ref2, 1e-300));
}

static void validate(const CylindricalModalResult& r, std::size_t count) {
    require(r.converged, "ellipticity convergence solve did not converge");
    require(r.frequencies_hz.size() == count, "wrong number of convergence frequencies");
    require(r.relative_residual < 2e-4, "ellipticity convergence residual too large");
    for (std::size_t i = 0; i < count; ++i) {
        require(std::isfinite(r.frequencies_hz[i]) && r.frequencies_hz[i] > 0.0,
                "invalid convergence frequency");
        if (i) require(r.frequencies_hz[i] >= r.frequencies_hz[i - 1],
                       "convergence frequencies are not sorted");
    }
}

int main() {
    constexpr double E = 1.0e11;
    constexpr double nu = 0.30;
    constexpr double rho = 7800.0;
    constexpr double q = 0.85;
    const std::vector<int> harmonics{0, 1, 2, 3};
    constexpr std::size_t modes = 4;

    // Meridional refinement: the same physical q and Fourier basis are used;
    // only the Q4 meridional mesh changes.
    const auto mesh2 = make_patch(2, 2);
    const auto mesh4 = make_patch(4, 4);
    const auto mesh8 = make_patch(8, 8);

    const auto solve = [&](const BellMeshData& mesh, unsigned angular_samples,
                           const std::vector<int>& hs) {
        const auto system = assemble_cylindrical_global_system(
            mesh, E, nu, rho, hs, angular_samples, 0.0, q);
        return solve_cylindrical_global_modes(
            system, clamp_lower_edge(mesh, system.dofs_per_node),
            static_cast<int>(modes), 300, 5e-5, true);
    };

    const auto r2 = solve(mesh2, 64, harmonics);
    const auto r4 = solve(mesh4, 64, harmonics);
    const auto r8 = solve(mesh8, 64, harmonics);
    validate(r2, modes); validate(r4, modes); validate(r8, modes);

    const double refine24 = max_relative_change(r2.frequencies_hz, r4.frequencies_hz);
    const double refine48 = max_relative_change(r4.frequencies_hz, r8.frequencies_hz);
    // These are convergence gates for the synthetic patch, not accuracy claims
    // for production bells. The refined solution must settle substantially.
    require(refine48 < refine24, "meridional refinement did not reduce spectral change");
    require(refine48 < 0.08, "8x8 refinement still changes low spectrum excessively");

    // Angular quadrature convergence on the finest meridional mesh.
    const auto system32 = assemble_cylindrical_global_system(
        mesh8, E, nu, rho, harmonics, 32, 0.0, q);
    const auto system64 = assemble_cylindrical_global_system(
        mesh8, E, nu, rho, harmonics, 64, 0.0, q);
    const auto system128 = assemble_cylindrical_global_system(
        mesh8, E, nu, rho, harmonics, 128, 0.0, q);
    const double matrix_k3264 = sparse_relative_difference(system32.stiffness, system64.stiffness);
    const double matrix_k64128 = sparse_relative_difference(system64.stiffness, system128.stiffness);
    const double matrix_m3264 = sparse_relative_difference(system32.mass, system64.mass);
    const double matrix_m64128 = sparse_relative_difference(system64.mass, system128.mass);
    require(matrix_k64128 < matrix_k3264, "angular stiffness quadrature did not settle");
    require(matrix_m64128 < matrix_m3264, "angular mass quadrature did not settle");
    require(matrix_k64128 < 1e-8 && matrix_m64128 < 1e-8,
            "128-sample angular quadrature matrix change remains excessive");

    const auto a32 = solve(mesh8, 32, harmonics);
    const auto a64 = solve(mesh8, 64, harmonics);
    const auto a128 = solve(mesh8, 128, harmonics);
    validate(a32, modes); validate(a64, modes); validate(a128, modes);
    const double ang3264 = max_relative_change(a32.frequencies_hz, a64.frequencies_hz);
    const double ang64128 = max_relative_change(a64.frequencies_hz, a128.frequencies_hz);
    // The matrices are the primary quadrature convergence contract. The
    // eigenvalue changes are additionally bounded, but their monotonicity is
    // not used as a hard gate because iterative eigen-solvers can introduce
    // platform-dependent roundoff at the 1e-10..1e-8 level.
    require(ang64128 < 2e-3, "128-sample angular quadrature changed low spectrum excessively");

    // Harmonic enrichment: adding Fourier content should settle toward the
    // coupled non-axisymmetric solution. Frequency rank is used only for this
    // convergence diagnostic; no modal-family identity claim is made here.
    const auto h01 = solve(mesh8, 64, {0, 1});
    const auto h012 = solve(mesh8, 64, {0, 1, 2});
    const auto h0123 = solve(mesh8, 64, {0, 1, 2, 3});
    validate(h01, modes); validate(h012, modes); validate(h0123, modes);
    const double harm12 = max_relative_change(h01.frequencies_hz, h012.frequencies_hz);
    const double harm23 = max_relative_change(h012.frequencies_hz, h0123.frequencies_hz);
    require(harm23 < harm12, "harmonic enrichment did not settle the low spectrum");
    require(harm23 < 0.03, "four-harmonic low spectrum remains too sensitive to enrichment");

    std::cout << "cylindrical_ellipticity_convergence: mesh/angular/harmonic refinement PASS\n"
              << "mesh_2x2_to_4x4_max_rel=" << refine24
              << " mesh_4x4_to_8x8_max_rel=" << refine48 << '\n'
              << "angular_32_to_64_max_rel=" << ang3264
              << " angular_64_to_128_max_rel=" << ang64128 << '\n'
              << "angular_matrix_K_32_to_64_rel=" << matrix_k3264
              << " angular_matrix_K_64_to_128_rel=" << matrix_k64128 << '\n'
              << "angular_matrix_M_32_to_64_rel=" << matrix_m3264
              << " angular_matrix_M_64_to_128_rel=" << matrix_m64128 << '\n'
              << "harmonics_01_to_012_max_rel=" << harm12
              << " harmonics_012_to_0123_max_rel=" << harm23 << '\n'
              << "q=" << q << " modes=" << modes << '\n';
    return 0;
}
