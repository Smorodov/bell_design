#include "fem/cylindrical_global_modes.hpp"
#include "fem/cylindrical_modal_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i)
        if (std::abs(mesh.nodes[i].z) < 1e-12)
            for (std::size_t d = 0; d < dpn; ++d)
                fixed.push_back(static_cast<int>(i * dpn + d));
    return fixed;
}

static void validate(const CylindricalModalResult& result, std::size_t modes) {
    require(result.converged, "modal family solve did not converge");
    require(result.frequencies_hz.size() >= modes && result.mode_shapes.size() >= modes,
            "modal family solve returned too few modes");
    require(result.relative_residual < 2e-4, "modal family residual too large");
}

int main() {
    constexpr double E = 1.0e11;
    constexpr double nu = 0.30;
    constexpr double rho = 7800.0;
    const auto mesh = make_patch(4, 4);
    const std::vector<int> harmonics{0, 1, 2, 3};
    constexpr std::size_t pair_dim = 2;

    // q=1 is the rotationally symmetric reference. The lowest pair is
    // degenerate, so its individual eigenvectors are not physically unique.
    const auto reference_system = assemble_cylindrical_global_system(
        mesh, E, nu, rho, harmonics, 64, 0.0, 1.0);
    const auto fixed = clamp_lower_edge(mesh, reference_system.dofs_per_node);
    const auto reference = solve_cylindrical_global_modes(
        reference_system, fixed, 6, 1800, 1e-4, true);
    validate(reference, 6);

    const double pair_split = std::abs(reference.frequencies_hz[1] - reference.frequencies_hz[0]) /
                              std::max(reference.frequencies_hz[0], 1e-300);
    require(pair_split < 1e-3, "q=1 reference pair is not sufficiently degenerate");

    // The reference mass is deliberately kept as a fixed common inner product.
    // This makes overlap a tracking coordinate, not a second q-dependent
    // physical quantity. It also makes the test independent of eigenvector
    // normalization and sign.
    const auto& metric_mass = reference_system.mass;
    const double self_subspace = mass_weighted_subspace_similarity(
        reference, 0, reference, 0, pair_dim, metric_mass);
    require(self_subspace > 1.0 - 1e-10, "reference degenerate subspace is not self-consistent");

    std::cout << "q=1 spectrum:"; for (std::size_t i=0;i<6;++i) std::cout << " " << reference.frequencies_hz[i]; std::cout << "\n";
    std::cout << "q=1 pair_Hz=" << reference.frequencies_hz[0] << ','
              << reference.frequencies_hz[1] << " relative_split=" << pair_split << '\n';

    double previous_subspace = 1.0;
    for (double q : {0.95, 0.85, 0.70}) {
        const auto system = assemble_cylindrical_global_system(
            mesh, E, nu, rho, harmonics, 64, 0.0, q);
        const auto candidate = solve_cylindrical_global_modes(
            system, fixed, 6, 1800, 1e-4, true);
        validate(candidate, 6);

        double best_subspace = -1.0;
        std::size_t best_i = 0, best_j = 1;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = i + 1; j < 6; ++j) {
                const double score = mass_weighted_subspace_similarity(
                    reference, 0, candidate, i, 1, metric_mass);
                const double score2 = mass_weighted_subspace_similarity(
                    reference, 1, candidate, j, 1, metric_mass);
                const double alt1 = mass_weighted_subspace_similarity(
                    reference, 0, candidate, j, 1, metric_mass);
                const double alt2 = mass_weighted_subspace_similarity(
                    reference, 1, candidate, i, 1, metric_mass);
                const double assignment = std::max(std::min(score, score2), std::min(alt1, alt2));
                if (assignment > best_subspace) {
                    best_subspace = assignment;
                    best_i = i; best_j = j;
                }
            }
        }
        const double ov00 = mass_weighted_modal_overlap(reference, 0, candidate, best_i, metric_mass);
        const double ov01 = mass_weighted_modal_overlap(reference, 0, candidate, best_j, metric_mass);
        const double ov10 = mass_weighted_modal_overlap(reference, 1, candidate, best_i, metric_mass);
        const double ov11 = mass_weighted_modal_overlap(reference, 1, candidate, best_j, metric_mass);
        const double pair_similarity = mass_weighted_subspace_similarity(
            reference, {0, 1}, candidate, {best_i, best_j}, metric_mass);
        const double direct_best = std::max(ov00 + ov11, ov01 + ov10) / 2.0;
        const double split = std::abs(candidate.frequencies_hz[1] - candidate.frequencies_hz[0]) /
                             std::max(candidate.frequencies_hz[0], 1e-300);

        // The physical pair must remain identifiable as a two-dimensional
        // family even when its two individual vectors rotate inside the pair.
        require(pair_similarity > 0.60, "ellipticity destroyed the reference modal-family identity");
        require(pair_similarity <= 1.0 + 1e-9, "subspace similarity exceeded unity");
        require(split > pair_split + 5e-4, "ellipticity did not split the reference pair");

        // As q moves away from the circular limit, the family should evolve
        // continuously. A mild decrease is expected; a sudden collapse would
        // indicate that sorted eigenvalue rank is no longer a safe family label.
        require(pair_similarity > previous_subspace - 0.20,
                "modal-family overlap collapsed abruptly during q sweep");
        previous_subspace = pair_similarity;

        std::cout << "q=" << q
                  << " pair_Hz=" << candidate.frequencies_hz[0] << ',' << candidate.frequencies_hz[1]
                  << " split=" << split
                  << " tracked_pair=" << best_i << "," << best_j
                  << " pair_similarity=" << pair_similarity
                  << " best_individual_overlap=" << direct_best
                  << " residual=" << candidate.relative_residual << '\n';
    }

    // Sign reversal of one candidate eigenvector must not change overlap.
    auto sign_flipped = reference;
    for (double& v : sign_flipped.mode_shapes[1]) v = -v;
    const double sign_invariant = mass_weighted_modal_overlap(
        reference, 1, sign_flipped, 1, metric_mass);
    require(sign_invariant > 1.0 - 1e-10, "modal overlap is not invariant to eigenvector sign");

    std::cout << "cylindrical_modal_family_tracking: degenerate-subspace and ellipticity sweep PASS\n";
    return 0;
}
