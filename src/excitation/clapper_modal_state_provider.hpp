#pragma once

#include "core/bell_types.h"
#include "core/modal_types.hpp"
#include "clapper_modal_evaluator.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"

namespace BellOptimizer::Excitation {

// Owns the complete FEM -> inner-wall contact -> released modal acoustic state
// pipeline for one candidate. All contact geometry is SI; BellGeometry remains
// in the optimizer's native units and is converted by the FEM model adapter.
struct ClapperModalStateProviderConfig {
    BellMeshConfig mesh;
    FEM::CylindricalModalEvaluationConfig fem;
    ClapperModalEvaluationConfig clapper;
};

struct ClapperModalState {
    ModalData modal_data;
    std::vector<Point2D> inner_profile_si;
    ClapperModalEvaluationResult excitation;
};

// Evaluates the clapper against modal data and the matching inner profile already
// produced by the candidate's final-geometry FEM solve. This is the production
// integration entry point: it performs no mesh generation and no FEM solve.
// The caller must ensure that modes and profile belong to the same final geometry
// and use SI coordinates for the profile.
bool evaluateClapperModalStateFromCandidate(
    const ClapperModalEvaluationConfig& config,
    const std::vector<Point2D>& inner_profile_si,
    const std::vector<ModalMode>& final_geometry_modes,
    ClapperModalEvaluationResult& result);

// Standalone diagnostic convenience path. It owns FEM model construction and
// solve, so it must not be called from an optimizer candidate evaluation.
// Returns false on invalid geometry, unconverged FEM solve, incomplete modal
// state, failed contact integration, or inadmissible strike.
bool evaluateClapperModalStateForGeometry(
    const BellGeometry& geometry,
    const BellMaterial& material,
    const ClapperModalStateProviderConfig& config,
    ClapperModalState& result);

} // namespace BellOptimizer::Excitation
