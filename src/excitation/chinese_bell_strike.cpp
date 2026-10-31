#include "chinese_bell_strike.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BellOptimizer::Excitation {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;

struct FourierWeight {
    std::vector<double> cos;
    std::vector<double> sin;
};

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

FourierWeight angular_projection(const std::vector<int>& harmonics, double center,
                                 double width, unsigned samples) {
    FourierWeight out;
    out.cos.assign(harmonics.size(), 0.0);
    out.sin.assign(harmonics.size(), 0.0);
    double total = 0.0;
    for (unsigned i = 0; i < samples; ++i) {
        const double phi = 2.0 * PI * static_cast<double>(i) / static_cast<double>(samples);
        const double d = wrap_angle(phi - center);
        const double w = std::exp(-0.5 * (d / width) * (d / width));
        total += w;
        for (std::size_t k = 0; k < harmonics.size(); ++k) {
            const int h = harmonics[k];
            if (h == 0) out.cos[k] += w;
            else {
                out.cos[k] += w * std::cos(static_cast<double>(h) * phi);
                out.sin[k] += w * std::sin(static_cast<double>(h) * phi);
            }
        }
    }
    if (!(total > 0.0) || !std::isfinite(total)) return {};
    for (std::size_t k = 0; k < harmonics.size(); ++k) {
        out.cos[k] /= total;
        out.sin[k] /= total;
    }
    return out;
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

    const auto fw = angular_projection(harmonics, config.azimuth_rad,
                                       config.angular_width_rad, angular_samples);
    if (fw.cos.size() != harmonics.size()) return {};

    std::vector<SurfaceFrame> frames(n);
    std::vector<double> nodal_area(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        frames[i] = surface_frame(mesh.outer_profile, i, ellipticity, config.azimuth_rad);

    // Integrate surface area segment-by-segment. Each segment contributes half
    // of its area to each endpoint (lumped trapezoidal surface quadrature).
    // This avoids the old centered-difference weighting, which over-weighted
    // interior profile nodes relative to the two endpoints.
    const double c = std::cos(config.azimuth_rad);
    const double s = std::sin(config.azimuth_rad);
    const double inv_q2 = 1.0 / (ellipticity * ellipticity);
    const double dphi_metric = c * c + s * s * inv_q2;
    const double alpha = 1.0 / std::sqrt(dphi_metric);
    const double log_alpha_derivative =
        -s * c * (inv_q2 - 1.0) / dphi_metric;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const auto& a = mesh.outer_profile[i];
        const auto& b = mesh.outer_profile[i + 1];
        const double dr = b.r - a.r;
        const double dz = b.z - a.z;
        const double r_mid = std::max(0.0, 0.5 * (a.r + b.r));
        const double metric = std::sqrt(
            dz * dz * (1.0 + log_alpha_derivative * log_alpha_derivative)
            + alpha * alpha * dr * dr);
        const double segment_area = alpha * r_mid * metric;
        nodal_area[i] += 0.5 * segment_area;
        nodal_area[i + 1] += 0.5 * segment_area;
    }

    double norm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = (mesh.outer_profile[i].z - config.center_z) / config.meridional_width;
        norm += std::exp(-0.5 * d * d) * nodal_area[i];
    }
    if (!(norm > std::numeric_limits<double>::epsilon()) || !std::isfinite(norm)) return {};

    for (std::size_t i = 0; i < n; ++i) {
        const double d = (mesh.outer_profile[i].z - config.center_z) / config.meridional_width;
        const double wz = std::exp(-0.5 * d * d) * nodal_area[i] / norm;
        const auto& frame = frames[i];
        const std::size_t node_base = i * 3u * nb;
        std::size_t b = 0;
        for (std::size_t k = 0; k < harmonics.size(); ++k) {
            const double scale_c = wz * fw.cos[k];
            out.dof_values[node_base + 3u * b + 0u] = frame.nr * scale_c;
            out.dof_values[node_base + 3u * b + 1u] = frame.nphi * scale_c;
            out.dof_values[node_base + 3u * b + 2u] = frame.nz * scale_c;
            ++b;
            if (harmonics[k] != 0) {
                const double scale_s = wz * fw.sin[k];
                out.dof_values[node_base + 3u * b + 0u] = frame.nr * scale_s;
                out.dof_values[node_base + 3u * b + 1u] = frame.nphi * scale_s;
                out.dof_values[node_base + 3u * b + 2u] = frame.nz * scale_s;
                ++b;
            }
        }
    }
    return out;
}

} // namespace BellOptimizer::Excitation
