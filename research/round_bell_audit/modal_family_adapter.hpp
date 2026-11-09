#pragma once

#include "../../src/core/modal_types.hpp"
#include <cstddef>
#include <vector>

namespace BellOptimizer::Research {

struct ModalFamilyAdapterConfig {
    // Adjacent frequencies within this relative gap are treated as one
    // numerically degenerate eigenspace. Keep deliberately tight by default.
    double degeneracy_relative_tolerance = 1.0e-5;
    double minimum_harmonic_share = 0.50;
};

struct ModalFamilyGroup {
    std::size_t first_mode = 0;
    std::size_t mode_count = 0;
    double representative_frequency_hz = 0.0;
    int dominant_azimuthal_order = -1;
    double dominant_harmonic_share = 0.0;
    double combined_strike_projection = 0.0;
    double projection_energy = 0.0;
    std::vector<double> harmonic_energy_shares; // m=0,1,... in basis order
};

// Computes mass-lumped harmonic energy shares and combines adjacent degenerate
// modes. Projection norm is invariant to orthogonal rotations/sign flips within
// an exactly degenerate, mass-orthonormal eigenspace. This is a research adapter;
// it intentionally does not assign western partial names or alter production.
bool buildModalFamilyGroups(const ModalData& modal_data,
                            const std::vector<double>& signed_strike_projections,
                            const ModalFamilyAdapterConfig& config,
                            std::vector<ModalFamilyGroup>& groups);

} // namespace BellOptimizer::Research
