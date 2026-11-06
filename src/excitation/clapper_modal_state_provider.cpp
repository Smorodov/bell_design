#include "clapper_modal_state_provider.hpp"
#include "fem/cylindrical_bell_model.hpp"

#include <cmath>
#include <exception>

namespace BellOptimizer::Excitation {

bool evaluateClapperModalStateFromCandidate(
    const ClapperModalEvaluationConfig& config,
    const std::vector<Point2D>& inner_profile_si,
    const std::vector<ModalMode>& final_geometry_modes,
    ClapperModalEvaluationResult& result) {
    result = {};
    if (inner_profile_si.size() < 2 || final_geometry_modes.size() < 5)
        return false;
    return evaluateClapperModalState(config, inner_profile_si, final_geometry_modes, result)
        && result.strike_height_admissible;
}

bool evaluateClapperModalStateForGeometry(
    const BellGeometry& geometry,
    const BellMaterial& material,
    const ClapperModalStateProviderConfig& config,
    ClapperModalState& result) {
    result = {};
    try {
        const auto model = FEM::make_cylindrical_bell_model(geometry, material, config.mesh);
        ModalInfluence neutral;
        neutral.radial.assign(model.mesh_si.nodes.size(), 0.0);
        neutral.axial.assign(model.mesh_si.nodes.size(), 0.0);
        if (!FEM::evaluate_cylindrical_bell_modes(model, neutral, config.fem, result.modal_data))
            return false;
        result.inner_profile_si = model.mesh_si.inner_profile;
        if (!evaluateClapperModalStateFromCandidate(
                config.clapper, result.inner_profile_si, result.modal_data.modes, result.excitation))
            return false;
        return true;
    } catch (const std::exception&) {
        result = {};
        return false;
    }
}

} // namespace BellOptimizer::Excitation
