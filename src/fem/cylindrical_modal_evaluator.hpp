#pragma once

#include "cylindrical_bell_model.hpp"
#include "core/modal_types.hpp"

namespace BellOptimizer::FEM {

struct CylindricalModalEvaluationConfig {
    std::vector<int> harmonics{0, 1, 2, 3};
    unsigned angular_samples = 64;
    int requested_modes = 10;
    int max_iterations = 1200;
    double tolerance = 1.0e-7;
    // The optimizer canonicalizes geometry to H=1000 mm before direct FEM.
    // This gate removes residual rigid-body numerical modes below the production
    // minimum frequency without hiding physical modes in the canonical search space.
    double minimum_frequency_hz = 1.0;
    bool use_rcm_ordering = true;
};

// Direct production-side modal evaluation for the vector cylindrical FEM path.
// It deliberately contains no coordinate morphing or frequency post-correction.
// The returned amplitude is a mass-projected meridional influence diagnostic.
// Experimental Chinese-bell strike logic lives under research/chinese_bell and
// is not part of the production library.
bool evaluate_cylindrical_bell_modes(
    const CylindricalBellModel& model,
    const ModalInfluence& influence,
    const CylindricalModalEvaluationConfig& config,
    ModalData& result);

} // namespace BellOptimizer::FEM
