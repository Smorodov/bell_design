#pragma once

#include "bell_types.h"
#include "modal_types.hpp"

#include <vector>

namespace BellOptimizer::Excitation {

struct LocalizedStrikeConfig {
    double center_z = 0.0;          // meridional coordinate, same units as mesh
    double width = 20.0;            // Gaussian 1-sigma width along the profile
    double normal_fraction = 1.0;   // 1 = outward normal, 0 = tangent
    double tangent_fraction = 0.0;
};

// Computes unit-total ring-strike modal impulse projections. The spatial
// impulse density is Gaussian along the outer profile and uniform in azimuth.
// The result is normalized so that the integral of |j| over the patch is one;
// therefore total impulse is an external scalar and is not hidden here.
std::vector<double> modalImpulse(const BellMeshData& mesh,
                                 const std::vector<ModalMode>& modes,
                                 const LocalizedStrikeConfig& config);

} // namespace BellOptimizer::Excitation
