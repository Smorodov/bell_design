#pragma once

#include "coordinate_morphing.hpp"
#include "bell_types.h"
#include "modal_types.hpp"

namespace BellOptimizer::Math {

struct MechanicalMorphingResult {
    double base_mass = 0.0;
    double base_strain_energy = 0.0;
    double morphed_mass = 0.0;
    double morphed_strain_energy = 0.0;
    double base_frequency_hz = 0.0;
    double morphed_frequency_hz = 0.0;
    double relative_mass_change = 0.0;
    double relative_energy_change = 0.0;
    double relative_frequency_change = 0.0;
};

// Evaluate the Rayleigh quotient of an axisymmetric FEM mode after an affine
// spatial morph x'=Sx. The original axisymmetric displacement is reconstructed
// in 3-D cylindrical coordinates, transformed with the same S, and integrated
// over azimuth. This is a reduced-order morphing calculation, not a new FEM
// solution for the morphed geometry.
MechanicalMorphingResult evaluate_mechanical_morphing(
    const BellMeshData& mesh,
    const BellMaterial& material,
    const ModalMode& mode,
    const CoordinateScale3D& scale);

} // namespace BellOptimizer::Math
