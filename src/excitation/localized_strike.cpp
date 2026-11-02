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

bool modalImpulseCircular(const BellMeshData& mesh,
                          const std::vector<ModalMode>& modes,
                          const CircularLocalizedStrikeConfig& config,
                          std::vector<double>& result) {
    result.clear();
    const auto& profile = mesh.outer_profile;
    if (profile.size() < 2 || modes.empty() ||
        !std::isfinite(config.center_z) ||
        !(config.meridional_sigma > 0.0) || !std::isfinite(config.meridional_sigma) ||
        !(config.circumferential_sigma > 0.0) || !std::isfinite(config.circumferential_sigma) ||
        !std::isfinite(config.azimuth_rad) ||
        !std::isfinite(config.normal_fraction) || !std::isfinite(config.tangent_fraction))
        return false;

    const double direction_norm = std::hypot(config.normal_fraction, config.tangent_fraction);
    if (!(direction_norm > 0.0) || !std::isfinite(direction_norm)) return false;
    const double dn = config.normal_fraction / direction_norm;
    const double dt = config.tangent_fraction / direction_norm;

    const std::size_t point_count = profile.size();
    const std::size_t basis_count = modes.front().outer_surface_harmonics.size();
    if (basis_count == 0 ||
        point_count > std::numeric_limits<std::size_t>::max() / basis_count) return false;
    const std::size_t expected = point_count * basis_count;

    for (const auto& mode : modes) {
        if (mode.outer_surface_harmonics.size() != basis_count ||
            mode.outer_surface_fourier_radial.size() != expected ||
            mode.outer_surface_fourier_azimuthal.size() != expected ||
            mode.outer_surface_fourier_axial.size() != expected) return false;
        bool has_zero = false;
        for (std::size_t b = 0; b < basis_count; ++b) {
            const int h = mode.outer_surface_harmonics[b];
            if (h == 0) has_zero = true;
            if (h == std::numeric_limits<int>::min()) return false;
            for (std::size_t prior = 0; prior < b; ++prior)
                if (mode.outer_surface_harmonics[prior] == h) return false;
            if (h != 0 && std::find(mode.outer_surface_harmonics.begin(),
                                   mode.outer_surface_harmonics.end(), -h) ==
                              mode.outer_surface_harmonics.end()) return false;
            for (std::size_t k = 0; k < point_count; ++k) {
                const std::size_t index = k * basis_count + b;
                if (!std::isfinite(mode.outer_surface_fourier_radial[index]) ||
                    !std::isfinite(mode.outer_surface_fourier_azimuthal[index]) ||
                    !std::isfinite(mode.outer_surface_fourier_axial[index])) return false;
            }
        }
        if (!has_zero) return false;
    }

    struct Segment {
        double r = 0.0;
        double z = 0.0;
        double ds = 0.0;
        double nr = 0.0;
        double nz = 0.0;
        double tr = 0.0;
        double tz = 0.0;
        double meridional_weight = 0.0;
        std::size_t profile_segment = 0;
    };
    std::vector<Segment> segments;
    segments.reserve(point_count - 1u);
    double normalization = 0.0;
    for (std::size_t k = 0; k + 1u < point_count; ++k) {
        const double dr = profile[k + 1u].r - profile[k].r;
        const double dz = profile[k + 1u].z - profile[k].z;
        const double ds = std::hypot(dr, dz);
        const double r = 0.5 * (profile[k].r + profile[k + 1u].r);
        const double z = 0.5 * (profile[k].z + profile[k + 1u].z);
        if (!(ds > 0.0) || !(r > 0.0) || !std::isfinite(ds) || !std::isfinite(r) || !std::isfinite(z))
            continue;
        const double tr = dr / ds;
        const double tz = dz / ds;
        const double nr = tz;
        const double nz = -tr;
        const double normalized_z = (z - config.center_z) / config.meridional_sigma;
        const double weight = std::exp(-0.5 * normalized_z * normalized_z);
        if (!(weight > 0.0) || !std::isfinite(weight)) continue;
        segments.push_back({r, z, ds, nr, nz, -tz, tr, weight, k});
        // The circumferential Gaussian is specified by physical arc length.
        // Its angular integral is proportional to sigma_s/r, which cancels
        // the surface Jacobian r in dS = r*ds*dphi. The common factor
        // sqrt(2*pi)*sigma_s cancels in the unit-impulse normalization.
        normalization += ds * weight;
    }
    if (!(normalization > std::numeric_limits<double>::epsilon()) || !std::isfinite(normalization))
        return false;

    result.assign(modes.size(), 0.0);
    for (std::size_t mode_index = 0; mode_index < modes.size(); ++mode_index) {
        const auto& mode = modes[mode_index];
        double projection = 0.0;
        for (const auto& x : segments) {
            const double dir_r = x.nr * dn + x.tr * dt;
            const double dir_z = x.nz * dn + x.tz * dt;
            const double area_weight = x.ds * x.meridional_weight / normalization;
            const double angular_sigma = config.circumferential_sigma / x.r;
            for (std::size_t b = 0; b < basis_count; ++b) {
                const int label = mode.outer_surface_harmonics[b];
                const int order = label < 0 ? -label : label;
                double angular_weight = 1.0;
                if (order > 0) {
                    const double m = static_cast<double>(order);
                    const double attenuation = std::exp(-0.5 * m * m * angular_sigma * angular_sigma);
                    const double phase = m * std::remainder(config.azimuth_rad, 2.0 * PI);
                    angular_weight = attenuation * (label > 0 ? std::cos(phase) : std::sin(phase));
                }
                const std::size_t i0 = x.profile_segment * basis_count + b;
                const std::size_t i1 = (x.profile_segment + 1u) * basis_count + b;
                const double ur = 0.5 * (mode.outer_surface_fourier_radial[i0] +
                                         mode.outer_surface_fourier_radial[i1]);
                const double uz = 0.5 * (mode.outer_surface_fourier_axial[i0] +
                                         mode.outer_surface_fourier_axial[i1]);
                projection += area_weight * angular_weight * (ur * dir_r + uz * dir_z);
            }
        }
        if (!std::isfinite(projection)) { result.clear(); return false; }
        result[mode_index] = projection;
    }
    return true;
}

} // namespace BellOptimizer::Excitation
