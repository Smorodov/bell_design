#include "internal_clapper_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {

bool find_spherical_tip_contact_geometry(
    const std::vector<Point2D>& profile,
    double centreR, double centreZ, double radius,
    SphericalTipContactGeometry& result) {
    result = {};
    if (profile.size() < 2 || !std::isfinite(centreR) || centreR < 0.0 ||
        !std::isfinite(centreZ) || !std::isfinite(radius) || radius <= 0.0)
        return false;
    for (const auto& p : profile)
        if (!std::isfinite(p.r) || p.r < 0.0 || !std::isfinite(p.z)) return false;

    const double rr = std::numeric_limits<double>::infinity();
    double bestD2 = rr;
    SphericalTipContactGeometry best;
    for (std::size_t i = 0; i + 1 < profile.size(); ++i) {
        const auto& a = profile[i];
        const auto& b = profile[i + 1];
        const double dr = b.r - a.r;
        const double dz = b.z - a.z;
        const double length2 = dr * dr + dz * dz;
        if (!(length2 > 1.0e-24) || !std::isfinite(length2)) continue;
        const double t = std::clamp(((centreR - a.r) * dr + (centreZ - a.z) * dz) / length2, 0.0, 1.0);
        const double pr = a.r + t * dr;
        const double pz = a.z + t * dz;
        const double nr = centreR - pr;
        const double nz = centreZ - pz;
        const double d2 = nr * nr + nz * nz;
        if (d2 < bestD2) {
            const double distance = std::sqrt(d2);
            if (!(distance > 1.0e-14)) continue; // normal is undefined at zero distance
            // Profile order is crown-to-mouth. For this orientation, (-dz, dr)
            // points from the metal wall into the bell cavity. Reject a sphere
            // centre lying on the exterior side of this local segment.
            const double tangentLength = std::sqrt(length2);
            const double cavityNr = -dz / tangentLength;
            const double cavityNz = dr / tangentLength;
            const double side = nr * cavityNr + nz * cavityNz;
            if (!(side > 0.0)) continue;
            bestD2 = d2;
            best.valid = true;
            best.segment_index = i;
            best.segment_parameter = t;
            best.surface_point = {pr, pz};
            best.outward_normal = {cavityNr, cavityNz};
            best.centre_distance = distance;
            best.gap = distance - radius;
        }
    }
    if (!best.valid) return false;
    result = best;
    return true;
}

} // namespace BellOptimizer::Excitation
