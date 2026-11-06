#include "core/bell_types.h"
#include "core/bell_material.h"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

struct HarmonicSignature {
    int order = -1;
    double dominant_fraction = 0.0;
};

static HarmonicSignature dominant_harmonic(const ModalMode& mode) {
    const std::size_t basis_count = mode.fem_harmonics.size();
    require(basis_count > 0, "FEM harmonic labels are missing");
    require(mode.fem_displacement.size() == mode.fem_mass_lumped_weights.size(),
            "FEM displacement/mass layouts differ");
    require(mode.fem_displacement.size() % (3u * basis_count) == 0,
            "FEM displacement layout is malformed");

    std::vector<double> energy(basis_count, 0.0);
    for (std::size_t i = 0; i < mode.fem_displacement.size(); ++i) {
        const std::size_t basis = i % basis_count;
        const double u = mode.fem_displacement[i];
        const double w = mode.fem_mass_lumped_weights[i];
        require(std::isfinite(u) && std::isfinite(w) && w >= 0.0,
                "non-finite FEM modal field or mass weight");
        energy[basis] += u * u * w;
    }
    double total = 0.0;
    for (double e : energy) total += e;
    require(total > 0.0 && std::isfinite(total), "modal harmonic energy is zero");
    const auto best = std::max_element(energy.begin(), energy.end());
    const std::size_t index = static_cast<std::size_t>(best - energy.begin());
    return {std::abs(mode.fem_harmonics[index]), *best / total};
}

int main() {
    BellMaterial material;
    BellMeshConfig mesh_config;
    mesh_config.nNodes = 15;
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    geometry.ellipticity = 1.0;

    const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh_config);
    FEM::CylindricalModalEvaluationConfig fem;
    fem.harmonics = {0, 1, 2, 3};
    fem.angular_samples = 64;
    fem.requested_modes = 32;
    fem.max_iterations = 5000;
    fem.tolerance = 1.0e-6;
    fem.minimum_frequency_hz = 1.0;
    Excitation::MeridionalInfluence influence;
    ModalData modal_data;
    require(FEM::evaluate_cylindrical_bell_modes(
                model, influence.make(model.mesh_si), fem, modal_data),
            "circular bell vector-FEM modal evaluation failed");
    require(modal_data.modes.size() >= 20, "too few modes for partial-family audit");

    std::cout << "index,frequency_hz,dominant_abs_harmonic,dominant_mass_fraction,amplitude\n";
    std::vector<HarmonicSignature> signatures;
    signatures.reserve(modal_data.modes.size());
    for (std::size_t i = 0; i < modal_data.modes.size(); ++i) {
        const auto signature = dominant_harmonic(modal_data.modes[i]);
        signatures.push_back(signature);
        std::cout << i << ',' << std::setprecision(12)
                  << modal_data.modes[i].frequency_hz << ',' << signature.order << ','
                  << signature.dominant_fraction << ',' << modal_data.modes[i].amplitude << '\n';
    }

    // This only detects candidate pairs from frequency and dominant harmonic
    // content. It deliberately reports mixed-harmonic eigenspaces instead of
    // treating a dominant label as a proven physical family identity.
    constexpr double pair_tolerance = 5.0e-5;
    std::size_t identified_pairs = 0;
    std::size_t ambiguous_pairs = 0;
    for (std::size_t i = 0; i + 1 < modal_data.modes.size(); ++i) {
        const auto& a = modal_data.modes[i];
        const auto& b = modal_data.modes[i + 1];
        const double gap = std::abs(a.frequency_hz - b.frequency_hz) /
            std::max(a.frequency_hz, b.frequency_hz);
        if (gap <= pair_tolerance && signatures[i].order > 0 &&
            signatures[i].order == signatures[i + 1].order) {
            ++identified_pairs;
            if (signatures[i].dominant_fraction <= 0.95 ||
                signatures[i + 1].dominant_fraction <= 0.95) ++ambiguous_pairs;
        }
    }
    require(identified_pairs >= 3, "too few near-degenerate modal pairs detected");
    require(ambiguous_pairs > 0, "expected mixed-harmonic eigenspaces were not detected");
    std::cout << "near_degenerate_pairs_with_matching_dominant_order=" << identified_pairs << '\n';
    std::cout << "pairs_with_ambiguous_harmonic_content=" << ambiguous_pairs << '\n';
    std::cout << "clapper_modal_spectrum_family_audit: PASS\n";
    return 0;
}
