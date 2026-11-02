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

// A physically localized patch on a circular/axisymmetric outer wall. Both
// widths are 1-sigma lengths in the same units as mesh.outer_profile. The
// circumferential Gaussian is periodized on the circle and projected
// analytically onto the retained Fourier basis; no angular quadrature is used.
// The result is a unit-total-impulse generalized-force projection per mode.
struct CircularLocalizedStrikeConfig {
    double center_z = 0.0;
    double meridional_sigma = 0.005;
    double circumferential_sigma = 0.005;
    double azimuth_rad = 0.0;
    double normal_fraction = 1.0;
    double tangent_fraction = 0.0;
};

// Returns false for incomplete/inconsistent Fourier boundary data or invalid
// geometry/configuration. On success, result has one signed projection per
// input mode. The mode's Fourier labels use +m for cos(m*phi), -m for
// sin(m*phi), and 0 for the axisymmetric component.
bool modalImpulseCircular(const BellMeshData& mesh,
                          const std::vector<ModalMode>& modes,
                          const CircularLocalizedStrikeConfig& config,
                          std::vector<double>& result);

} // namespace BellOptimizer::Excitation
