#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "fem/cylindrical_modal_field.hpp"
#include "core/modal_tracking.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

struct SolvedBell {
    CylindricalBellModel model;
    std::vector<int> harmonics;
    ModalData modes;
};

static BellGeometry make_geometry() {
    BellGeometry g;
    g.H = 658.3175709034034;
    g.R_mouth = 355.9758748202059;
    g.t_top = 22.85675347179014;
    g.t_bottom = 46.62524667849419;
    g.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    g.ellipticity = 0.85;
    return g;
}

static SolvedBell solve(int nodes, std::vector<int> harmonics, unsigned angular_samples) {
    BellMaterial material;
    BellMeshConfig mesh;
    mesh.nNodes = nodes;
    auto geometry = make_geometry();
    geometry.scale(1000.0 / geometry.H);

    SolvedBell result;
    result.model = make_cylindrical_bell_model(geometry, material, mesh);
    result.harmonics = std::move(harmonics);

    CylindricalModalEvaluationConfig cfg;
    cfg.harmonics = result.harmonics;
    cfg.angular_samples = angular_samples;
    cfg.requested_modes = 14;
    cfg.max_iterations = 5000;
    cfg.tolerance = 2.0e-6;
    cfg.minimum_frequency_hz = 1.0;

    Excitation::MeridionalInfluence influence;
    require(evaluate_cylindrical_bell_modes(result.model,
                influence.make(result.model.mesh_si), cfg, result.modes),
            "direct FEM solve failed");
    require(result.modes.modes.size() >= 12, "too few modes returned for continuation matrix");
    for (const auto& mode : result.modes.modes) {
        require(mode.fem_displacement.size() == result.model.mesh_si.nodes.size() *
                    mode.fem_harmonics.size() * 3u,
                "invalid full FEM field export");
    }
    return result;
}

struct OwnedField {
    std::vector<int> harmonics;
    const BellMeshData* mesh = nullptr;
    const std::vector<double>* shape = nullptr;
    double ellipticity = 1.0;
    double density = 0.0;

    CylindricalModalField view() const {
        return {mesh, &harmonics, shape, ellipticity, density};
    }
};

static OwnedField field(const SolvedBell& s, std::size_t mode_index) {
    const auto& m = s.modes.modes.at(mode_index);
    OwnedField out;
    // The modal-field API accepts nonnegative harmonic orders; the solver's
    // exported labels alternate cosine (+h) and sine (-h) basis entries.
    for (int h : m.fem_harmonics) if (h >= 0) out.harmonics.push_back(h);
    out.mesh = &s.model.mesh_si;
    out.shape = &m.fem_displacement;
    out.ellipticity = s.model.ellipticity;
    out.density = s.model.density_kg_m3;
    return out;
}

struct PairMatch {
    double similarity = 0.0;
    std::size_t a0 = 0, a1 = 0, b0 = 0, b1 = 0;
    double common_coverage = 0.0;
};

static PairMatch best_pair_match(const SolvedBell& a, const SolvedBell& b,
                                 std::size_t candidate_begin = 8,
                                 std::size_t candidate_end = 14) {
    PairMatch best;
    const std::size_t na = a.modes.modes.size();
    const std::size_t nb = std::min(candidate_end, b.modes.modes.size());
    // The Chinese-bell readiness study identified the 11th/12th modes
    // (zero-based 10/11) as the nearly-degenerate reference family. Keep that
    // family fixed and search nearby candidate branches after discretization changes.
    if (na <= 11 || nb <= candidate_begin) return best;
    const std::size_t i = 10, j = 11;
    {
        auto fa0 = field(a, i);
        auto fa1 = field(a, j);
        const std::vector<CylindricalModalField> fa{fa0.view(), fa1.view()};
        for (std::size_t k = candidate_begin; k < nb; ++k)
            for (std::size_t l = k + 1; l < nb; ++l) {
            auto fb0 = field(b, k);
            auto fb1 = field(b, l);
            const std::vector<CylindricalModalField> fb{fb0.view(), fb1.view()};
            double coverage = 0.0;
            const double score = cylindrical_modal_subspace_similarity(fa, fb, 32, 0.98, &coverage);
            if (score > best.similarity) best = {score, i, j, k, l, coverage};
        }
    }
    return best;
}

static void print_match(const char* label, const PairMatch& m) {
    std::cout << label << " similarity=" << m.similarity
              << " coverage=" << m.common_coverage
              << " ref_pair=" << m.a0 << ',' << m.a1
              << " candidate_pair=" << m.b0 << ',' << m.b1 << '\n';
}

int main() {
    try {
        const std::vector<int> base{0, 2, 4};
        const std::vector<int> enriched{0, 2, 4, 6};

        // Step 1: identical geometry and Fourier basis, different profile meshes.
        const auto mesh_ref = solve(11, base, 32);
        const auto mesh_candidate = solve(15, base, 32);
        const auto mesh_match = best_pair_match(mesh_ref, mesh_candidate);
        print_match("mesh_only", mesh_match);
        require(mesh_match.similarity > 0.70,
                "physical modal subspace did not continue under mesh refinement");
        require(mesh_match.common_coverage > 0.98,
                "mesh refinement did not preserve common physical domain");

        // Step 2: fixed mesh and geometry, enriched Fourier basis.
        const auto fourier_candidate = solve(11, enriched, 32);
        const auto fourier_match = best_pair_match(mesh_ref, fourier_candidate);
        print_match("fourier_only", fourier_match);
        require(fourier_match.similarity > 0.50,
                "no physically compatible pair found after Fourier enrichment");
        require(fourier_match.common_coverage > 0.98,
                "Fourier enrichment unexpectedly changed the physical domain");

        // Step 3: mesh and Fourier basis change simultaneously.
        const auto combined_candidate = solve(15, enriched, 32);
        const auto combined_match = best_pair_match(mesh_ref, combined_candidate);
        print_match("mesh_and_fourier", combined_match);
        require(combined_match.similarity > 0.45,
                "no physically compatible pair found under simultaneous mesh/Fourier change");
        require(combined_match.common_coverage > 0.98,
                "combined refinement unexpectedly changed the physical domain");

        // Compare the existing boundary matcher as a baseline on the mesh-only case.
        auto baseline_candidate = mesh_candidate.modes;
        const auto boundary_baseline = trackModalFamilies(mesh_ref.modes, baseline_candidate);
        std::cout << "boundary_baseline used_shape=" << (boundary_baseline.used_shape_tracking ? 1 : 0)
                  << " min_similarity=" << boundary_baseline.minimum_similarity << '\n';

        std::cout << "cylindrical_modal_field_fem_matrix: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "cylindrical_modal_field_fem_matrix: FAIL: " << e.what() << '\n';
        return 1;
    }
}
