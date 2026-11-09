#pragma once

#include "modal_family_adapter.hpp"
#include "../../src/excitation/localized_strike.hpp"

namespace BellOptimizer::Research {

// Single hand-off object for the optimizer integration boundary. It contains
// research diagnostics only; it is deliberately not an objective-function term.
struct RoundBellCandidateAnalysis {
    std::vector<double> signed_strike_projections;
    std::vector<ModalFamilyGroup> modal_families;
    double maximum_absolute_strike_projection = 0.0;
};

// Projects a circular localized strike onto the supplied direct-FEM modes and
// groups the resulting modes. The mesh, modal fields and mode ordering must be
// from the same FEM solve. On failure, result is cleared and false is returned.
bool analyzeRoundBellCandidate(const BellMeshData& mesh,
                               const ModalData& modal_data,
                               const Excitation::CircularLocalizedStrikeConfig& strike,
                               const ModalFamilyAdapterConfig& adapter_config,
                               RoundBellCandidateAnalysis& result);

} // namespace BellOptimizer::Research
