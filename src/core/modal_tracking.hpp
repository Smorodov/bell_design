#pragma once

#include "modal_types.hpp"

#include <vector>

namespace BellOptimizer {

struct ModalTrackingReport {
    bool used_shape_tracking = false;
    bool ambiguous = false;
    double minimum_similarity = 0.0;
    double mean_similarity = 0.0;
    std::size_t matched_modes = 0;
};

// Reorders candidate modes so that modal-family identity follows normalized
// outer-wall displacement shape rather than frequency rank alone. The reference
// and candidate must have compatible outer-profile sampling. If shape data are
// unavailable, the candidate is left in frequency order and the report says so.
// Returns the minimum principal-angle cosine between two boundary-trace
// subspaces. Each basis vector must use the same flattened radial/axial
// sampling. Sign, scale and orthogonal rotations within either basis do not
// affect the result. Returns -1 when the inputs are invalid or rank deficient.
double boundaryTraceSubspaceSimilarity(
    const std::vector<std::vector<double>>& reference_basis,
    const std::vector<std::vector<double>>& candidate_basis);

ModalTrackingReport trackModalFamilies(const ModalData& reference,
                                       ModalData& candidate,
                                       double ambiguity_similarity = 0.05,
                                       double near_degenerate_frequency_ratio = 0.01);

} // namespace BellOptimizer
