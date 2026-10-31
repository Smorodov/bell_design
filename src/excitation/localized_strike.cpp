#include "localized_strike.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;

struct Sample {
    double z = 0.0;
    double r = 0.0;
    double ds = 0.0;
    double nr = 0.0;
    double nz = 0.0;
    double tr = 0.0;
    double tz = 0.0;
    double w = 0.0;
};

std::vector<Sample> samples(const BellMeshData& mesh, const LocalizedStrikeConfig& c) {
    const auto& p = mesh.outer_profile;
    std::vector<Sample> out;
    if (p.size() < 2 || !(c.width > 0.0) || !std::isfinite(c.width)) return out;

    const double nf = c.normal_fraction;
    const double tf = c.tangent_fraction;
    const double norm = std::hypot(nf, tf);
    if (!(norm > 0.0) || !std::isfinite(norm)) return out;
    // Normalize the direction independently from the scalar impulse density.
    // Otherwise (normal_fraction, tangent_fraction) would silently change the
    // total impulse magnitude, even though the caller requested a unit-total
    // impulse.
    for (size_t i = 0; i + 1 < p.size(); ++i) {
        const double dr = p[i + 1].r - p[i].r;
        const double dz = p[i + 1].z - p[i].z;
        const double ds = std::hypot(dr, dz);
        if (!(ds > 0.0)) continue;
        const double r = 0.5 * (p[i].r + p[i + 1].r);
        const double z = 0.5 * (p[i].z + p[i + 1].z);
        // Profile is traversed from crown toward mouth. The outward normal
        // is the right-hand normal in the (r,z) plane.
        const double tr = dr / ds;
        const double tz = dz / ds;
        const double nr = tz;
        const double nz = -tr;
        const double shape = std::exp(-0.5 * std::pow((z - c.center_z) / c.width, 2.0));
        out.push_back({z, r, ds, nr, nz, -tz, tr, shape});
    }
    return out;
}

} // namespace

std::vector<double> modalImpulse(const BellMeshData& mesh,
                                 const std::vector<ModalMode>& modes,
                                 const LocalizedStrikeConfig& config) {
    std::vector<double> result(modes.size(), 0.0);
    const auto s = samples(mesh, config);
    if (s.empty()) return result;
    const double direction_norm = std::hypot(config.normal_fraction, config.tangent_fraction);
    if (!(direction_norm > 0.0) || !std::isfinite(direction_norm)) return result;
    const double dn = config.normal_fraction / direction_norm;
    const double dt = config.tangent_fraction / direction_norm;

    // Surface element of an axisymmetric ring: dS = 2*pi*r*ds.
    // Normalize the scalar spatial impulse density to unit total impulse.
    double total = 0.0;
    for (const auto& x : s) total += 2.0 * PI * std::max(0.0, x.r) * x.ds * x.w;
    if (!(total > std::numeric_limits<double>::epsilon()) || !std::isfinite(total)) return result;

    for (size_t m = 0; m < modes.size(); ++m) {
        const auto& mode = modes[m];
        if (mode.outer_surface_radial.size() != mesh.outer_profile.size() ||
            mode.outer_surface_axial.size() != mesh.outer_profile.size()) {
            continue;
        }
        // Interpolate the modal displacement linearly at each segment midpoint.
        for (size_t i = 0; i + 1 < mesh.outer_profile.size(); ++i) {
            const auto& x = s[i];
            const double ur = 0.5 * (mode.outer_surface_radial[i] + mode.outer_surface_radial[i + 1]);
            const double uz = 0.5 * (mode.outer_surface_axial[i] + mode.outer_surface_axial[i + 1]);
            const double dir_r = x.w * (x.nr * dn + x.tr * dt);
            const double dir_z = x.w * (x.nz * dn + x.tz * dt);
            result[m] += (ur * dir_r + uz * dir_z) * (2.0 * PI * std::max(0.0, x.r) * x.ds) / total;
        }
    }
    return result;
}

} // namespace BellOptimizer::Excitation
