#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"
#include "math_constants.hpp"

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

// Exact axisymmetric surface integral for a linear displacement field on one
// profile segment. z_ref is the acoustic origin.
void segment_moments(const Point2D& p0, const Point2D& p1,
                     double ur0, double uz0, double ur1, double uz1,
                     double z_ref, double& qd, double& dz_moment) {
    const double dr = p1.r - p0.r;
    const double dz = p1.z - p0.z;
    const double a0 = p0.r * (ur0 * dz - uz0 * dr);
    const double a1 = dr * (ur0 * dz - uz0 * dr)
                    + p0.r * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
    const double a2 = dr * ((ur1 - ur0) * dz - (uz1 - uz0) * dr);
    const double i0 = a0 + 0.5 * a1 + a2 / 3.0;
    const double i1 = (p0.z - z_ref) * i0
                    + dz * (0.5 * a0 + a1 / 3.0 + a2 / 4.0);
    qd += 2.0 * bell::kPi * i0;
    dz_moment += 2.0 * bell::kPi * i1;
}

bool synthetic_checks(const std::vector<Point2D>& profile) {
    if (profile.size() < 2) return false;
    const double zref = profile.back().z;

    double q = 0.0, d = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i)
        segment_moments(profile[i], profile[i + 1], 1.0, 0.0, 1.0, 0.0,
                        zref, q, d);

    // A uniform radial field has a nonzero monopole and a finite first moment.
    if (!(std::isfinite(q) && std::isfinite(d) && std::abs(q) > 1e-12)) return false;

    // Construct a continuous linear radial field ur = z + a and choose a so
    // that its net monopole is zero. A nonzero first moment must remain. This
    // is the cleanest synthetic check that the dipole term is not merely a
    // rescaled monopole.
    double qz = 0.0, dz0 = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i)
        segment_moments(profile[i], profile[i + 1], profile[i].z, 0.0,
                        profile[i + 1].z, 0.0, zref, qz, dz0);
    const double a = -qz / q;
    double qzero = 0.0, dzero = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i)
        segment_moments(profile[i], profile[i + 1], profile[i].z + a, 0.0,
                        profile[i + 1].z + a, 0.0, zref, qzero, dzero);
    const double zero_monopole_rel = std::abs(qzero) / std::max(1.0, std::abs(qz));
    std::cout << "zero_monopole_relative_residual=" << zero_monopole_rel
              << " zero_monopole_dipole=" << dzero << '\n';
    if (!(zero_monopole_rel < 1e-12 && std::abs(dzero) > 1e-10)) return false;

    // Scaling the displacement field must scale both moments by exactly the
    // same factor.
    double q2scale = 0.0, d2scale = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i)
        segment_moments(profile[i], profile[i + 1], 2.0, 0.0, 2.0, 0.0,
                        zref, q2scale, d2scale);
    const double scale_err = std::max(
        std::abs(q2scale - 2.0 * q) / std::max(1.0, std::abs(q)),
        std::abs(d2scale - 2.0 * d) / std::max(1.0, std::abs(d)));
    std::cout << "moment_scaling_relative_error=" << scale_err << '\n';
    if (scale_err > 1e-12) return false;

    // Translating the origin changes D by -delta*z*Q. This is an exact
    // multipole consistency relation, and catches sign/orientation errors.
    const double delta = 17.0;
    double q2 = 0.0, d2 = 0.0;
    for (size_t i = 0; i + 1 < profile.size(); ++i)
        segment_moments(profile[i], profile[i + 1], 1.0, 0.0, 1.0, 0.0,
                        zref + delta, q2, d2);
    const double expected = d - delta * q;
    const double rel = std::abs(d2 - expected) / std::max(1.0, std::abs(expected));
    std::cout << "origin_translation_relative_error=" << rel << '\n';
    return rel < 1e-12;
}

} // namespace

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(reference_geometry());
    if (!generator.buildMesh()) return 2;
    if (!synthetic_checks(generator.getOuterProfile())) return 3;

    BellOptimizer::InternalModalSolver solver(true, 10);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), material, ".", "multipole_audit", &modal_influence_1};
    if (!solver.solve(problem, data) || data.modes.size() < 6) return 4;

    constexpr double c = 343.0;
    // The compact-source criterion is based on the source geometry, not on
    // the potentially ill-conditioned ratio Dz/Qd. Coordinates are in mm.
    double z_min = generator.getOuterProfile().front().z;
    double z_max = z_min;
    double r_max = 0.0;
    for (const auto& point : generator.getOuterProfile()) {
        z_min = std::min(z_min, point.z);
        z_max = std::max(z_max, point.z);
        r_max = std::max(r_max, std::abs(point.r));
    }
    const double axial_span_m = (z_max - z_min) * 1e-3;
    const double diameter_m = 2.0 * r_max * 1e-3;
    const double characteristic_length_m = std::max(axial_span_m, diameter_m);
    if (!(std::isfinite(characteristic_length_m) && characteristic_length_m > 0.0)) return 9;

    std::cout << std::setprecision(12);
    std::cout << "geometry_axial_span_m=" << axial_span_m << '\n';
    std::cout << "geometry_max_diameter_m=" << diameter_m << '\n';
    std::cout << "geometry_characteristic_length_m=" << characteristic_length_m << '\n';
    std::cout << "mode,frequency_hz,Qd,Dz_moment_mm,abs_Qd_over_maxQd,Dz_over_Qd_mm,Dz_over_Qd_m,first_order_ratio_axial,first_order_ratio_45deg,kL_max_dimension,source_regime,monopole_status\n";

    size_t weak_monopoles = 0;
    size_t undefined_ratios = 0;
    size_t finite_count = 0;
    double max_q = 0.0;
    for (const auto& m : data.modes) max_q = std::max(max_q, std::abs(m.acoustic_volume_displacement));
    if (!(std::isfinite(max_q) && max_q > 1e-12)) return 6;

    // Ratios are ill-conditioned when Qd is near zero. Do not silently report
    // zero in that case: mark the ratio undefined and preserve Dz for review.
    const double ratio_floor = std::max(1e-14, max_q * 1e-8);
    for (const auto& m : data.modes) {
        const double q = m.acoustic_volume_displacement;
        const double d = m.acoustic_dipole_moment_z;
        if (!std::isfinite(q) || !std::isfinite(d) || !std::isfinite(m.frequency_hz)) return 5;
        ++finite_count;
        const double q_relative = std::abs(q) / max_q;
        const bool weak = q_relative < 0.02;
        if (weak) ++weak_monopoles;
        const double k = 2.0 * bell::kPi * m.frequency_hz / c;
        const double kL = k * characteristic_length_m;
        // Conservative regimes: kL <= 0.3 is compact; <= 1 is transitional.
        // Above 1, a low-order source multipole truncation is not justified
        // from compactness alone. This is not a far-field-distance test.
        const char* regime = kL <= 0.3 ? "COMPACT" :
                             (kL <= 1.0 ? "TRANSITIONAL" : "NONCOMPACT" );
        std::cout << m.mode_number << ',' << m.frequency_hz << ',' << q << ',' << d
                  << ',' << q_relative << ',';
        if (std::abs(q) > ratio_floor) {
            // BellGeometry and FEM coordinates are in millimetres. Convert the
            // first-moment length to metres before multiplying by the wave
            // number (1/m). This ratio is diagnostic, not a synthesis gain.
            const double length_mm = std::abs(d / q);
            const double length_m = length_mm * 1e-3;
            const double ratio_axial = k * length_m;
            const double ratio_45 = ratio_axial * std::sqrt(0.5);
            if (!std::isfinite(length_mm) || !std::isfinite(length_m) ||
                !std::isfinite(ratio_axial)) return 7;
            std::cout << length_mm << ',' << length_m << ',' << ratio_axial << ','
                      << ratio_45 << ',' << kL << ',' << regime << ','
                      << (weak ? "WEAK_MONOPOLE" : "MONOPOLE_NOT_SUPPRESSED") << '\n';
        } else {
            ++undefined_ratios;
            std::cout << "NA,NA,NA,NA," << kL << "," << regime << ",RATIO_ILL_CONDITIONED\n";
        }
    }

    std::cout << "finite_modes=" << finite_count << '\n';
    std::cout << "weak_monopole_count=" << weak_monopoles << '\n';
    std::cout << "ill_conditioned_ratio_count=" << undefined_ratios << '\n';
    std::cout << "reference_sound_speed_mps=" << c << '\n';
    std::cout << "compact_source_count=" << std::count_if(data.modes.begin(), data.modes.end(), [&](const auto& m) { return (2.0 * bell::kPi * m.frequency_hz / c) * characteristic_length_m <= 0.3; }) << '\n';
    std::cout << "noncompact_source_count=" << std::count_if(data.modes.begin(), data.modes.end(), [&](const auto& m) { return (2.0 * bell::kPi * m.frequency_hz / c) * characteristic_length_m > 1.0; }) << '\n';
    std::cout << "interpretation: kL uses the maximum of axial span and outer diameter; kL assesses source compactness, not listener far-field distance; first-order ratios are diagnostics only\n";
    return finite_count == data.modes.size() ? 0 : 8;
}
