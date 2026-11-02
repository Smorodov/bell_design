#pragma once

#include "core/bell_types.h"
#include <cstddef>
#include <vector>

namespace BellOptimizer::Excitation {

// Closest meridional contact between a spherical tip and an axisymmetric
// inner-wall profile. Coordinates and radius use the same unit (production
// FEM profiles are SI; BellGenerator profiles are millimetres).
struct SphericalTipContactGeometry {
    bool valid = false;
    std::size_t segment_index = 0;
    double segment_parameter = 0.0;
    Point2D surface_point{};
    Point2D outward_normal{}; // from bell material into the cavity, toward tip centre
    double centre_distance = 0.0;
    double gap = 0.0; // positive = separated, zero = tangent, negative = overlap
};

// Finds the closest point on the piecewise-linear inner meridian to the sphere
// centre. The result is a geometric gap only; no force law or dynamics is applied.
bool find_spherical_tip_contact_geometry(
    const std::vector<Point2D>& inner_profile,
    double tip_centre_r, double tip_centre_z, double tip_radius,
    SphericalTipContactGeometry& result);

} // namespace BellOptimizer::Excitation
