#include "chinese_bell_strike.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BellOptimizer::Excitation {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;

std::size_t basis_count(const std::vector<int>& harmonics) {
    if (harmonics.empty() || harmonics.front() != 0)
        throw std::invalid_argument("Chinese strike: harmonics must start at zero");
    std::size_t n = 0;
    int previous = -1;
    for (int h : harmonics) {
        if (h < 0 || h <= previous)
            throw std::invalid_argument("Chinese strike: harmonics must be increasing and unique");
        n += (h == 0 ? 1u : 2u);
        previous = h;
    }
    return n;
}

double wrap_angle(double a) {
    while (a > PI) a -= 2.0 * PI;
    while (a < -PI) a += 2.0 * PI;
    return a;
}

struct SurfaceFrame {
    double nr = 0.0;
    double nphi = 0.0;
    double nz = 0.0;
};

SurfaceFrame surface_frame(const std::vector<Point2D>& p, std::size_t i,
                           double ellipticity, double phi) {
    const std::size_t j0 = i == 0 ? 0 : i - 1;
    const std::size_t j1 = std::min(i + 1, p.size() - 1);
    const double dr = p[j1].r - p[j0].r;
    const double dz = p[j1].z - p[j0].z;
    const double c = std::cos(phi), s = std::sin(phi);
    const double inv_q2 = 1.0 / (ellipticity * ellipticity);
    const double d = c * c + s * s * inv_q2;
    const double alpha = 1.0 / std::sqrt(d);
    const double log_alpha_derivative = -s * c * (inv_q2 - 1.0) / d;

    // For r(phi,z)=alpha(phi)*r_ref(z), the outward normal in cylindrical
    // components is proportional to (dz, -dz*d(log(alpha))/d(phi), -alpha*dr).
    // This reduces exactly to the meridional normal when q=1.
    double nr = dz;
    double nphi = -dz * log_alpha_derivative;
    double nz = -alpha * dr;
    const double length = std::sqrt(nr * nr + nphi * nphi + nz * nz);
    if (length > 0.0) { nr /= length; nphi /= length; nz /= length; }
    else { nr = nphi = nz = 0.0; }

    return {nr, nphi, nz};
}

} // namespace

FEMStrikeInfluence make_chinese_bell_strike_influence(
    const BellMeshData& mesh,
    const std::vector<int>& harmonics,
    double ellipticity,
    const ChineseBellStrikeConfig& config,
    unsigned angular_samples) {
    FEMStrikeInfluence out;
    if (mesh.nodes.empty() || mesh.outer_profile.size() < 2) return out;
    if (!(ellipticity > 0.0 && ellipticity <= 1.0) || !std::isfinite(ellipticity))
        throw std::invalid_argument("Chinese strike: ellipticity must be finite in (0,1]");
    if (!(config.meridional_width > 0.0) || !(config.angular_width_rad > 0.0) ||
        !std::isfinite(config.meridional_width) || !std::isfinite(config.angular_width_rad)) return out;
    if (angular_samples < 32 || angular_samples % 2 != 0) return out;

    const std::size_t nb = basis_count(harmonics);
    const std::size_t n = mesh.outer_profile.size();
    if (mesh.nodes.size() < n) return out;
    out.dof_values.assign(mesh.nodes.size() * 3u * nb, 0.0);
    out.center_z = config.center_z;
    out.azimuth_rad = config.azimuth_rad;
    out.id = "chinese_bell_localized_surface_normal_strike";

    // The generalized force coefficient for a Fourier displacement basis is
    // the actual surface integral
    //
    //   F_b = integral_S g(z,phi) n(phi) B_b(phi) dS .
    //
    // For an ellipse both the surface normal and dS vary with phi.  The old
    // implementation evaluated the normal and meridional surface weight only
    // at the strike azimuth and then multiplied by a scalar Fourier projection.
    // That is not the same integral and becomes physically inconsistent as
    // ellipticity departs from one.  Integrate the complete vector load over
    // phi instead; the same quadrature is then used by the FEM basis.
    std::vector<double> z_weight(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double d = (mesh.outer_profile[i].z - config.center_z) / config.meridional_width;
        z_weight[i] = std::exp(-0.5 * d * d);
    }

    auto nodal_area_per_phi = [&](std::size_t i, double phi) {
        if (i >= n) return 0.0;
        const double c = std::cos(phi), s = std::sin(phi);
        const double inv_q2 = 1.0 / (ellipticity * ellipticity);
        const double metric_phi = c * c + s * s * inv_q2;
        const double alpha = 1.0 / std::sqrt(metric_phi);
        const double log_alpha_derivative =
            -s * c * (inv_q2 - 1.0) / metric_phi;
        double area = 0.0;
        if (i > 0) {
            const auto& a = mesh.outer_profile[i - 1];
            const auto& b = mesh.outer_profile[i];
            const double dr = b.r - a.r;
            const double dz = b.z - a.z;
            const double r_mid = std::max(0.0, 0.5 * (a.r + b.r));
            const double metric = std::sqrt(
                dz * dz * (1.0 + log_alpha_derivative * log_alpha_derivative) +
                alpha * alpha * dr * dr);
            area += 0.5 * alpha * r_mid * metric;
        }
        if (i + 1 < n) {
            const auto& a = mesh.outer_profile[i];
            const auto& b = mesh.outer_profile[i + 1];
            const double dr = b.r - a.r;
            const double dz = b.z - a.z;
            const double r_mid = std::max(0.0, 0.5 * (a.r + b.r));
            const double metric = std::sqrt(
                dz * dz * (1.0 + log_alpha_derivative * log_alpha_derivative) +
                alpha * alpha * dr * dr);
            area += 0.5 * alpha * r_mid * metric;
        }
        return area;
    };

    const double dphi = 2.0 * PI / static_cast<double>(angular_samples);
    std::vector<double> force(out.dof_values.size(), 0.0);
    double total_force_weight = 0.0;
    for (unsigned ip = 0; ip < angular_samples; ++ip) {
        const double phi = dphi * static_cast<double>(ip);
        const double dphi_angle = wrap_angle(phi - config.azimuth_rad);
        const double angular_weight = std::exp(
            -0.5 * (dphi_angle / config.angular_width_rad) *
            (dphi_angle / config.angular_width_rad));
        // The normal is evaluated at the same profile node and azimuth as the
        // corresponding force quadrature point.
        
        std::vector<double> basis_values;
        basis_values.reserve(nb);
        for (int h : harmonics) {
            basis_values.push_back(h == 0 ? 1.0 : std::cos(static_cast<double>(h) * phi));
            if (h != 0) basis_values.push_back(std::sin(static_cast<double>(h) * phi));
        }
        for (std::size_t i = 0; i < n; ++i) {
            const double area = nodal_area_per_phi(i, phi);
            const auto frame = surface_frame(mesh.outer_profile, i, ellipticity, phi);
            const double scalar = z_weight[i] * angular_weight * area * dphi;
            if (!(scalar >= 0.0) || !std::isfinite(scalar)) return {};
            total_force_weight += scalar;
            const std::size_t base = i * 3u * nb;
            for (std::size_t b = 0; b < nb; ++b) {
                force[base + 3u * b + 0u] += scalar * frame.nr * basis_values[b];
                force[base + 3u * b + 1u] += scalar * frame.nphi * basis_values[b];
                force[base + 3u * b + 2u] += scalar * frame.nz * basis_values[b];
            }
        }
    }
    if (!(total_force_weight > std::numeric_limits<double>::epsilon()) ||
        !std::isfinite(total_force_weight)) return {};
    for (double& value : force) value /= total_force_weight;
    out.dof_values = std::move(force);
    return out;
}

} // namespace BellOptimizer::Excitation
