#include "fem/sparse_eigensolver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using BellOptimizer::detail::SparseMatrix;
using BellOptimizer::detail::SparseTriplet;
using BellOptimizer::detail::solve_lowest_sparse_modes;

static SparseMatrix make_neumann_chain(int n, double epsilon) {
    std::vector<SparseTriplet> entries;
    entries.reserve(static_cast<std::size_t>(3 * n));
    for (int i = 0; i < n; ++i) {
        const double diagonal = (i == 0 || i == n - 1) ? 1.0 : 2.0;
        entries.push_back({i, i, diagonal + epsilon});
        if (i > 0) entries.push_back({i, i - 1, -1.0});
        if (i + 1 < n) entries.push_back({i, i + 1, -1.0});
    }
    return SparseMatrix::from_triplets(n, std::move(entries));
}

static bool ic0_contract() {
    using BellOptimizer::detail::sparse_eigen_detail::IC0Preconditioner;
    constexpr int n = 24;
    std::vector<SparseTriplet> entries;
    for (int i = 0; i < n; ++i) {
        entries.push_back({i, i, 2.25});
        if (i > 0) entries.push_back({i, i - 1, -1.0});
        if (i + 1 < n) entries.push_back({i, i + 1, -1.0});
    }
    const auto a = SparseMatrix::from_triplets(n, std::move(entries));
    IC0Preconditioner preconditioner;
    if (!preconditioner.build(a, 0.0)) return false;
    std::vector<double> x_exact(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) x_exact[static_cast<std::size_t>(i)] = 1.0 + 0.25 * i;
    const auto b = a.multiply(x_exact);
    std::vector<double> x;
    if (!preconditioner.solve(b, x)) return false;
    double max_error = 0.0;
    for (int i = 0; i < n; ++i)
        max_error = std::max(max_error, std::abs(x[static_cast<std::size_t>(i)] - x_exact[static_cast<std::size_t>(i)]));
    std::cout << "ic0_contract_max_error=" << max_error << '\n';
    return max_error < 1e-12;
}

int main() {
    if (!ic0_contract()) {
        std::cerr << "IC(0) preconditioner contract failed\n";
        return 10;
    }
    constexpr int n = 200;
    constexpr double epsilon = 1e-12;
    std::vector<SparseTriplet> mass_entries;
    mass_entries.reserve(n);
    for (int i = 0; i < n; ++i) mass_entries.push_back({i, i, 1.0});

    const auto stiffness = make_neumann_chain(n, epsilon);
    const auto mass = SparseMatrix::from_triplets(n, std::move(mass_entries));

    // The unshifted chain has a condition number of order 1/epsilon. The test
    // deliberately exercises the difficult near-null mode. The solver must not
    // report convergence merely because Ritz eigenvalues have stabilized.
    const auto result = solve_lowest_sparse_modes(stiffness, mass, 5, 120, 1e-8);

    if (result.converged) {
        if (!(result.final_relative_residual < 1e-7) ||
            !std::isfinite(result.final_relative_residual) ||
            !std::isfinite(result.final_eigenvalue_change)) {
            std::cerr << "solver reported convergence with invalid diagnostics\n";
            return 1;
        }
    } else {
        if (!std::isfinite(result.final_relative_residual) ||
            !std::isfinite(result.final_eigenvalue_change) ||
            result.iterations != 120) {
            std::cerr << "non-converged result did not preserve finite diagnostics\n";
            return 2;
        }
    }

    // The smallest mode is a near-null constant vector. The normwise backward
    // error must remain meaningful even though both K*x and lambda*M*x are tiny.
    // The first eigenvalue is close to epsilon and the corresponding vector is
    // checked independently against the known null-space direction.
    if (result.eigenvalues.empty() ||
        !std::isfinite(result.eigenvalues.front()) ||
        std::abs(result.eigenvalues.front() - epsilon) > 5e-7) {
        std::cerr << "near-null eigenvalue diagnostic failed: "
                  << (result.eigenvalues.empty() ? -1.0 : result.eigenvalues.front()) << '\n';
        return 3;
    }

    if (result.converged) {
        const auto& v = result.vectors.front();
        const double mean = [&]() { double s = 0.0; for (double x : v) s += x; return s / v.size(); }();
        double max_deviation = 0.0;
        for (double x : v) max_deviation = std::max(max_deviation, std::abs(x - mean));
        if (!(max_deviation < 1e-7)) {
            std::cerr << "near-null vector is not sufficiently close to the constant mode: " << max_deviation << "\n";
            return 4;
        }
    }

    std::cout << "sparse_conditioning_audit=PASS\n"
              << "n=" << n << " epsilon=" << epsilon
              << " converged=" << result.converged
              << " iterations=" << result.iterations
              << " first_eigenvalue=" << result.eigenvalues.front()
              << " final_residual=" << result.final_relative_residual
              << " final_eigenvalue_change=" << result.final_eigenvalue_change << '\n';
    return 0;
}
