#pragma once

#include "listener_geometry.hpp"
#include "../core/bell_types.h"
#include "../core/modal_types.hpp"

#include <complex>
#include <vector>

namespace BellOptimizer::Acoustics {

// Direct finite-distance evaluation of the existing outer-shell surface
// radiation proxy at one ideal point receiver. The result uses the same
// far-field-normalized units as the directional modal integral.
// This is a free-space Green-kernel approximation, not a complete exterior
// boundary-element solution for an unbaffled, open bell shell.
std::complex<double> evaluateFiniteDistanceRadiation(
    const ModalMode& mode,
    const std::vector<Point2D>& meridionalProfile,
    const ReceiverPoint& receiver,
    unsigned azimuthSamples = 192,
    double soundSpeedMPerS = 343.0,
    unsigned meridionalSubdivisions = 1);

} // namespace BellOptimizer::Acoustics
