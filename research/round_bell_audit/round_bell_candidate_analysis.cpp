#include "round_bell_candidate_analysis.hpp"
#include "../../src/excitation/localized_strike.hpp"

#include <algorithm>
#include <cmath>

namespace BellOptimizer::Research {

bool analyzeRoundBellCandidate(
    const BellMeshData& mesh,
    const ModalData& modal_data,
    const Excitation::CircularLocalizedStrikeConfig& strike,
    const ModalFamilyAdapterConfig& adapter_config,
    RoundBellCandidateAnalysis& result) {
    result = {};
    if (modal_data.modes.empty()) return false;

    if (!Excitation::modalImpulseCircular(
            mesh, modal_data.modes, strike, result.signed_strike_projections)) {
        result = {};
        return false;
    }

    for (double value : result.signed_strike_projections) {
        if (!std::isfinite(value)) {
            result = {};
            return false;
        }
        result.maximum_absolute_strike_projection =
            std::max(result.maximum_absolute_strike_projection, std::abs(value));
    }
    if (!(result.maximum_absolute_strike_projection > 1.0e-14) ||
        !buildModalFamilyGroups(modal_data, result.signed_strike_projections,
                                adapter_config, result.modal_families)) {
        result = {};
        return false;
    }
    return true;
}

} // namespace BellOptimizer::Research
