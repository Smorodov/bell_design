#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

BellGeometry reference_geometry() {
    BellGeometry g;
    g.H = 658.3175709034034;
    g.R_mouth = 355.9758748202059;
    g.t_top = 22.85675347179014;
    g.t_bottom = 46.62524667849419;
    g.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    return g;
}

// Exact segment flux for a linear displacement field on a linear profile.
double segment_qd(const Point2D& p0, const Point2D& p1,
                  double ur0, double uz0, double ur1, double uz1) {
    const double dr = p1.r - p0.r;
    const double dz = p1.z - p0.z;
    const double a0 = p0.r * (ur0 * dz - uz0 * dr);
    const double a1 = dr * (ur0 * dz - uz0 * dr)
                    + p0.r * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
    const double a2 = dr * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
    constexpr double PI = 3.1415926535897932384626433832795;
    return 2.0 * PI * (a0 + 0.5 * a1 + a2 / 3.0);
}

bool audit_surface_flux(const std::vector<Point2D>& profile) {
    if (profile.size() < 2) return false;

    // A uniform radial displacement has Qd = 2*pi*integral r dz for the
    // stored piecewise-linear profile. It is a useful orientation check.
    double q_radial = 0.0;
    double exact_radial = 0.0;
    // A uniform axial displacement has Qd = -pi*R_mouth^2 exactly.
    double q_axial = 0.0;
    // A linear radial displacement provides a nontrivial quadratic integrand.
    double q_linear = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i) {
        const auto& p0 = profile[i];
        const auto& p1 = profile[i + 1];
        q_radial += segment_qd(p0, p1, 1.0, 0.0, 1.0, 0.0);
        exact_radial += 2.0 * 3.1415926535897932384626433832795
                      * 0.5 * (p0.r + p1.r) * (p1.z - p0.z);
        q_axial  += segment_qd(p0, p1, 0.0, 1.0, 0.0, 1.0);
        q_linear += segment_qd(p0, p1, p0.r, 0.0, p1.r, 0.0);
    }

    const double R = profile.back().r;
    constexpr double PI = 3.1415926535897932384626433832795;
    const double exact_axial = -PI * R * R;
    const double radial_err = std::abs(q_radial - exact_radial) / std::max(1.0, std::abs(exact_radial));
    const double axial_err = std::abs(q_axial - exact_axial) / std::max(1.0, std::abs(exact_axial));

    // The same exact segment formula must be invariant to subdivision of a
    // segment. This isolates quadrature error from FE discretization error.
    const auto& p0 = profile.front();
    const auto& p1 = profile.back();
    const double base = segment_qd(p0, p1, p0.r, 0.25 * p0.r, p1.r, 0.25 * p1.r);
    double subdivided = 0.0;
    for (int j = 0; j < 16; ++j) {
        const double t0 = static_cast<double>(j) / 16.0;
        const double t1 = static_cast<double>(j + 1) / 16.0;
        Point2D a{p0.r + t0 * (p1.r - p0.r), p0.z + t0 * (p1.z - p0.z)};
        Point2D b{p0.r + t1 * (p1.r - p0.r), p0.z + t1 * (p1.z - p0.z)};
        subdivided += segment_qd(a, b, a.r, 0.25 * a.r, b.r, 0.25 * b.r);
    }
    const double subdivision_err = std::abs(base - subdivided) / std::max(1.0, std::abs(base));

    std::cout << "uniform_radial_relative_error=" << radial_err << '\n';
    std::cout << "q_radial=" << q_radial << " exact_radial=" << exact_radial << "\n";
    std::cout << "uniform_axial_relative_error=" << axial_err << '\n';
    std::cout << "subdivision_relative_error=" << subdivision_err << '\n';
    std::cout << "linear_field_qd=" << q_linear << '\n';
    return radial_err < 1e-12 && axial_err < 1e-12 && subdivision_err < 1e-12;
}

} // namespace

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(reference_geometry());
    if (!generator.buildMesh()) return 2;

    BellOptimizer::InternalModalSolver solver(true, 10);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "acoustic_audit", &modal_influence_1};
    if (!solver.solve(problem, data) || data.modes.size() < 6) return 3;

    if (!audit_surface_flux(generator.getOuterProfile())) return 4;

    size_t near_zero = 0;
    double max_abs_qd = 0.0;
    for (const auto& mode : data.modes) {
        if (!std::isfinite(mode.acoustic_volume_displacement) ||
            !std::isfinite(mode.modal_excitation)) return 5;
        max_abs_qd = std::max(max_abs_qd, std::abs(mode.acoustic_volume_displacement));
        if (std::abs(mode.acoustic_volume_displacement) < 0.02 * max_abs_qd && max_abs_qd > 0.0) ++near_zero;
    }

    std::cout << std::setprecision(12);
    std::cout << "modes=" << data.modes.size() << '\n';
    std::cout << "max_abs_modal_qd=" << max_abs_qd << '\n';
    std::cout << "near_zero_modal_qd_count=" << near_zero << '\n';

    // Small Qd is reported, not rejected: such modes are candidates for a
    // future multipole audit, not a reason to inject an empirical correction.
    return max_abs_qd > 1e-8 ? 0 : 6;
}
