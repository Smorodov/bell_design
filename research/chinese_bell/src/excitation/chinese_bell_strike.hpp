#pragma once

#include "core/bell_types.h"
#include <vector>

namespace BellOptimizer::Excitation {

struct ChineseBellStrikeConfig {
    double center_z = 0.0;          // meridional strike height, metres (SI)
    double azimuth_rad = 0.0;       // strike azimuth; 0=long axis, pi/2=short axis
    double meridional_width = 0.020; // Gaussian standard deviation along profile, metres (SI)
    double angular_width_rad = 0.18;// Gaussian width around azimuth
};

struct FEMStrikeInfluence {
    std::vector<double> dof_values; // full (node, component, Fourier basis) influence
    std::string id;
    double center_z = 0.0;
    double azimuth_rad = 0.0;
};

// Builds a localized unit-normal-load influence on the actual outer surface.
// The 3-D normal includes the azimuthal component introduced by ellipticity.
// Fourier coefficients are projected into the same cosine/sine basis as FEM.
FEMStrikeInfluence make_chinese_bell_strike_influence(
    const BellMeshData& mesh,
    const std::vector<int>& harmonics,
    double ellipticity,
    const ChineseBellStrikeConfig& config,
    unsigned angular_samples = 128);

} // namespace BellOptimizer::Excitation
