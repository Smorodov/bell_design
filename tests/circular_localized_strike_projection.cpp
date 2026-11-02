#include "core/modal_types.hpp"
#include "bell_generator.h"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_global_assembly.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "excitation/localized_strike.hpp"
#include "excitation/meridional_influence.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

BellOptimizer::ModalMode make_mode(const std::vector<int>& harmonics,
                                   double radial_cos, double radial_sin,
                                   double radial_zero) {
    BellOptimizer::ModalMode mode;
    mode.outer_surface_harmonics = harmonics;
    const std::size_t count = 2u * harmonics.size();
    mode.outer_surface_fourier_radial.assign(count, 0.0);
    mode.outer_surface_fourier_azimuthal.assign(count, 0.0);
    mode.outer_surface_fourier_axial.assign(count, 0.0);
    for (std::size_t b = 0; b < harmonics.size(); ++b) {
        double value = 0.0;
        if (harmonics[b] == 0) value = radial_zero;
        else if (harmonics[b] == 1) value = radial_cos;
        else if (harmonics[b] == -1) value = radial_sin;
        mode.outer_surface_fourier_radial[b] = value;
        mode.outer_surface_fourier_radial[harmonics.size() + b] = value;
    }
    return mode;
}
}

double evaluate_surface_component(const ModalMode& mode, std::size_t point,
                                  double azimuth, int component) {
    const std::size_t basis_count = mode.outer_surface_harmonics.size();
    const auto& values = component == 0 ? mode.outer_surface_fourier_radial
                       : component == 1 ? mode.outer_surface_fourier_azimuthal
                                        : mode.outer_surface_fourier_axial;
    double value = 0.0;
    for (std::size_t b = 0; b < basis_count; ++b) {
        const int harmonic = mode.outer_surface_harmonics[b];
        const double angular = harmonic == 0 ? 1.0
            : harmonic > 0 ? std::cos(static_cast<double>(harmonic) * azimuth)
                           : std::sin(static_cast<double>(-harmonic) * azimuth);
        value += values[point * basis_count + b] * angular;
    }
    return value;
}

int main() {
    try {
        BellMeshData mesh;
        // A unit-radius vertical meridian segment: the outward normal is +r.
        mesh.outer_profile = {{1.0, 0.0}, {1.0, 1.0}};

        const std::vector<int> basis{0, 1, -1};
        const auto cosine_mode = make_mode(basis, 1.0, 0.0, 0.0);
        const auto sine_mode = make_mode(basis, 0.0, 1.0, 0.0);
        const auto axisymmetric_mode = make_mode(basis, 0.0, 0.0, 1.0);
        const std::vector<ModalMode> modes{cosine_mode, sine_mode, axisymmetric_mode};

        CircularLocalizedStrikeConfig strike;
        strike.center_z = 0.5;
        strike.meridional_sigma = 0.2;
        strike.circumferential_sigma = 0.1;
        strike.azimuth_rad = 0.0;

        std::vector<double> at_zero;
        require(modalImpulseCircular(mesh, modes, strike, at_zero),
                "valid Fourier strike projection rejected");
        require(at_zero.size() == 3, "wrong projection count");

        const double attenuation = std::exp(-0.5 * 0.1 * 0.1);
        require(std::abs(at_zero[0] - attenuation) < 1e-12,
                "cosine mode projection at zero azimuth is incorrect");
        require(std::abs(at_zero[1]) < 1e-12,
                "sine mode projection at zero azimuth is not zero");
        require(std::abs(at_zero[2] - 1.0) < 1e-12,
                "axisymmetric projection should be azimuth invariant");

        const double theta = 0.73;
        strike.azimuth_rad = theta;
        std::vector<double> rotated;
        require(modalImpulseCircular(mesh, modes, strike, rotated),
                "rotated Fourier strike projection rejected");
        require(std::abs(rotated[0] - attenuation * std::cos(theta)) < 1e-12,
                "cosine modal projection does not obey rotation law");
        require(std::abs(rotated[1] - attenuation * std::sin(theta)) < 1e-12,
                "sine modal projection does not obey rotation law");
        require(std::abs(rotated[2] - at_zero[2]) < 1e-12,
                "m=0 projection changed under strike rotation");
        require(std::abs(std::hypot(rotated[0], rotated[1]) -
                         std::hypot(at_zero[0], at_zero[1])) < 1e-12,
                "degenerate cosine/sine pair norm changed under rotation");

        strike.circumferential_sigma = 0.2;
        std::vector<double> broader;
        require(modalImpulseCircular(mesh, modes, strike, broader),
                "broader contact projection rejected");
        require(std::abs(broader[0]) < std::abs(rotated[0]),
                "broader physical contact did not attenuate m=1");

        // The same analytic projection contract applies to higher modal
        // families: an isolated m=2 cosine/sine pair rotates by 2*theta and
        // preserves its pair norm. This is independent of eigenvector basis.
        ModalMode m2_cos, m2_sin;
        for (ModalMode* mode : {&m2_cos, &m2_sin}) {
            mode->outer_surface_harmonics = {0, 2, -2};
            mode->outer_surface_fourier_radial.assign(6, 0.0);
            mode->outer_surface_fourier_azimuthal.assign(6, 0.0);
            mode->outer_surface_fourier_axial.assign(6, 0.0);
        }
        m2_cos.outer_surface_fourier_radial[1] = 1.0;
        m2_cos.outer_surface_fourier_radial[4] = 1.0;
        m2_sin.outer_surface_fourier_radial[2] = 1.0;
        m2_sin.outer_surface_fourier_radial[5] = 1.0;
        CircularLocalizedStrikeConfig m2_strike = strike;
        m2_strike.azimuth_rad = 0.0;
        m2_strike.circumferential_sigma = 0.1;
        std::vector<double> m2_zero, m2_rotated;
        require(modalImpulseCircular(mesh, {m2_cos, m2_sin}, m2_strike, m2_zero),
                "m=2 zero-angle projection rejected");
        m2_strike.azimuth_rad = theta;
        require(modalImpulseCircular(mesh, {m2_cos, m2_sin}, m2_strike, m2_rotated),
                "m=2 rotated projection rejected");
        const double m2_attenuation = std::exp(-0.5 * 4.0 * 0.1 * 0.1);
        require(std::abs(m2_zero[0] - m2_attenuation) < 1e-12 &&
                std::abs(m2_zero[1]) < 1e-12,
                "m=2 basis projection is incorrect");
        require(std::abs(m2_rotated[0] - m2_attenuation * std::cos(2.0 * theta)) < 1e-12 &&
                std::abs(m2_rotated[1] - m2_attenuation * std::sin(2.0 * theta)) < 1e-12,
                "m=2 pair does not obey its doubled rotation phase");
        require(std::abs(std::hypot(m2_rotated[0], m2_rotated[1]) -
                         std::hypot(m2_zero[0], m2_zero[1])) < 1e-12,
                "m=2 pair norm changed under rotation");

        // A fixed physical arc-width footprint has area proportional to
        // meridional ds after integrating around phi: r cancels against
        // sigma_phi = sigma_s/r. Verify the projection does not reweight
        // equal-length segments by their circumference.
        BellMeshData varying_radius_mesh;
        varying_radius_mesh.outer_profile = {{1.0, 0.0}, {1.0, 1.0},
                                              {2.0, 1.0}, {2.0, 2.0}};
        ModalMode radial_patch_mode;
        radial_patch_mode.outer_surface_harmonics = {0};
        radial_patch_mode.outer_surface_fourier_radial = {1.0, 1.0, 0.0, 0.0};
        radial_patch_mode.outer_surface_fourier_azimuthal = {0.0, 0.0, 0.0, 0.0};
        radial_patch_mode.outer_surface_fourier_axial = {0.0, 0.0, 0.0, 0.0};
        CircularLocalizedStrikeConfig area_strike;
        area_strike.center_z = 1.0;
        area_strike.meridional_sigma = 1.0e6;
        area_strike.circumferential_sigma = 0.1;
        std::vector<double> area_projection;
        require(modalImpulseCircular(varying_radius_mesh, {radial_patch_mode},
                                     area_strike, area_projection),
                "variable-radius physical-footprint projection rejected");
        require(std::abs(area_projection[0] - 1.0 / 3.0) < 1.0e-10,
                "physical arc-width normalization incorrectly weights by radius");

        auto incomplete = modes;
        incomplete[0].outer_surface_harmonics = {0, 1};
        incomplete[0].outer_surface_fourier_radial.resize(4);
        incomplete[0].outer_surface_fourier_azimuthal.resize(4);
        incomplete[0].outer_surface_fourier_axial.resize(4);
        std::vector<double> rejected;
        require(!modalImpulseCircular(mesh, incomplete, strike, rejected),
                "unpaired Fourier basis must be rejected");

        // Integration and covariance audit: direct vector FEM supplies the
        // complete Fourier boundary trace. Check the exported layout against
        // the assembled [node][component][basis] global DOF convention.
        BellMaterial material;
        BellGeometry geometry;
        geometry.H = 658.3175709034034;
        geometry.R_mouth = 355.9758748202059;
        geometry.t_top = 22.85675347179014;
        geometry.t_bottom = 46.62524667849419;
        geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
        geometry.ellipticity = 1.0;
        geometry.scale(1000.0 / geometry.H);

        std::size_t finest_mode_count = 0;
        double coarsest_pair_frequency = 0.0;
        double finest_pair_frequency = 0.0;
        double coarsest_pair_power = 0.0;
        double finest_pair_power = 0.0;
        double maximum_rotation_error = 0.0;
        for (int profile_nodes : {15, 23, 31}) {
            BellMeshConfig mesh_config;
            mesh_config.nNodes = profile_nodes;
            const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh_config);
            FEM::CylindricalModalEvaluationConfig fem_config;
            fem_config.harmonics = {0, 1, 2};
            fem_config.angular_samples = 32;
            fem_config.requested_modes = 16;
            fem_config.max_iterations = 5000;
            fem_config.tolerance = 1.0e-6;
            fem_config.minimum_frequency_hz = 1.0;
            ModalData direct;
            MeridionalInfluence influence;
            require(FEM::evaluate_cylindrical_bell_modes(
                        model, influence.make(model.mesh_si), fem_config, direct),
                    "direct vector FEM failed on circular geometry");
            require(direct.modes.size() >= 2, "circular vector FEM returned too few modes");

            const auto& first = direct.modes[0];
            const auto& second = direct.modes[1];
            const double pair_split = std::abs(second.frequency_hz - first.frequency_hz) /
                                      std::max(first.frequency_hz, 1.0e-12);
            require(pair_split < 1.0e-6,
                    "lowest circular m=1 cosine/sine pair is not degenerate");
            const std::size_t basis_count = first.outer_surface_harmonics.size();
            require(basis_count == 5, "unexpected direct-FEM Fourier basis count");
            require(first.fem_displacement.size() ==
                        model.mesh_si.nodes.size() * 3u * basis_count,
                    "full FEM displacement has incorrect size");
            for (std::size_t k = 0; k < model.mesh_si.outer_profile.size(); ++k) {
                const std::size_t node_base = k * 3u * basis_count;
                require(std::abs(first.nodal_radial[k] -
                                 first.fem_displacement[node_base]) < 1.0e-12,
                        "n=0 radial nodal trace uses incorrect DOF index");
                require(std::abs(first.nodal_axial[k] -
                                 first.fem_displacement[node_base + 2u * basis_count]) < 1.0e-12,
                        "n=0 axial nodal trace uses incorrect DOF index");
                for (std::size_t b = 0; b < basis_count; ++b) {
                    const std::size_t dst = k * basis_count + b;
                    require(std::abs(first.outer_surface_fourier_radial[dst] -
                                     first.fem_displacement[node_base + b]) < 1.0e-12,
                            "radial Fourier trace does not match global FEM DOF layout");
                    require(std::abs(first.outer_surface_fourier_azimuthal[dst] -
                                     first.fem_displacement[node_base + basis_count + b]) < 1.0e-12,
                            "azimuthal Fourier trace does not match global FEM DOF layout");
                    require(std::abs(first.outer_surface_fourier_axial[dst] -
                                     first.fem_displacement[node_base + 2u * basis_count + b]) < 1.0e-12,
                            "axial Fourier trace does not match global FEM DOF layout");
                }
            }

            // Independently reconstruct the evaluator's generalized influence
            // using the assembled mass matrix and the documented component-major
            // Fourier layout. This catches an incorrect radial/axial DOF offset.
            const auto assembled = FEM::assemble_cylindrical_bell_system(
                model, fem_config.harmonics, fem_config.angular_samples);
            const auto mass_times_mode = assembled.mass.multiply(first.fem_displacement);
            const auto generalized_influence = influence.make(model.mesh_si);
            double expected_participation = 0.0;
            for (std::size_t node = 0; node < model.mesh_si.nodes.size(); ++node) {
                const std::size_t base = node * assembled.dofs_per_node;
                expected_participation += generalized_influence.radial[node] *
                    mass_times_mode[base] + generalized_influence.axial[node] *
                    mass_times_mode[base + 2u * basis_count];
            }
            require(std::abs(expected_participation - first.part_rz) <
                        1.0e-10 * std::max(1.0, std::abs(expected_participation)),
                    "generalized influence projection uses incorrect Fourier DOF offsets");

            CircularLocalizedStrikeConfig direct_strike;
            direct_strike.center_z = model.mesh_si.outer_profile.front().z +
                0.70 * (model.mesh_si.outer_profile.back().z -
                        model.mesh_si.outer_profile.front().z);
            direct_strike.meridional_sigma = 0.012;
            direct_strike.circumferential_sigma = 0.008;
            direct_strike.azimuth_rad = 0.0;
            std::vector<double> projection_zero;
            require(modalImpulseCircular(model.mesh_si, direct.modes, direct_strike,
                                         projection_zero),
                    "circular FEM zero-azimuth strike projection rejected");
            direct_strike.azimuth_rad = 0.63;
            std::vector<double> projection_rotated;
            require(modalImpulseCircular(model.mesh_si, direct.modes, direct_strike,
                                         projection_rotated),
                    "circular FEM rotated strike projection rejected");

            // Reconstruct the outer-wall modal velocity field from the complete
            // returned mode window. A rotated physical strike must rotate this
            // field, not merely change a scalar amplitude or pan balance.
            double error2 = 0.0;
            double reference2 = 0.0;
            constexpr int azimuth_samples = 48;
            const double theta = direct_strike.azimuth_rad;
            for (std::size_t k = 0; k < model.mesh_si.outer_profile.size(); ++k) {
                for (int ip = 0; ip < azimuth_samples; ++ip) {
                    const double phi = 2.0 * 3.14159265358979323846 * ip /
                                       static_cast<double>(azimuth_samples);
                    for (int component : {0, 2}) {
                        double rotated_field = 0.0;
                        double reference_field = 0.0;
                        for (std::size_t i = 0; i < direct.modes.size(); ++i) {
                            rotated_field += projection_rotated[i] *
                                evaluate_surface_component(direct.modes[i], k, phi, component);
                            reference_field += projection_zero[i] *
                                evaluate_surface_component(direct.modes[i], k, phi - theta, component);
                        }
                        error2 += (rotated_field - reference_field) *
                                  (rotated_field - reference_field);
                        reference2 += reference_field * reference_field;
                    }
                }
            }
            const double rotation_error = std::sqrt(error2 / std::max(reference2, 1.0e-30));
            maximum_rotation_error = std::max(maximum_rotation_error, rotation_error);
            require(rotation_error < 5.0e-4,
                    "FEM modal response violates circular strike rotation covariance");

            // The lowest m=1 eigenspace is a degenerate pair. Its summed
            // generalized-force power is independent of arbitrary cos/sin
            // eigenvector orientation and is recorded across profile meshes.
            const double pair_power = projection_zero[0] * projection_zero[0] +
                                      projection_zero[1] * projection_zero[1];
            require(std::isfinite(pair_power) && pair_power > 0.0,
                    "lowest circular modal pair has no strike response");
            if (profile_nodes == 15) {
                coarsest_pair_frequency = first.frequency_hz;
                coarsest_pair_power = pair_power;
            }
            if (profile_nodes == 31) {
                finest_pair_frequency = first.frequency_hz;
                finest_pair_power = pair_power;
            }
            finest_mode_count = direct.modes.size();
            std::cout << "mesh_nodes=" << profile_nodes
                      << " lowest_pair_frequency_Hz=" << first.frequency_hz
                      << " pair_power=" << pair_power
                      << " rotation_covariance_error=" << rotation_error << "\n";
        }
        require(coarsest_pair_frequency > 0.0 && finest_pair_frequency > 0.0 &&
                coarsest_pair_power > 0.0 && finest_pair_power > 0.0,
                "mesh audit did not collect its endpoint diagnostics");
        std::cout << "lowest_pair_frequency_change="
                  << (finest_pair_frequency - coarsest_pair_frequency) /
                     coarsest_pair_frequency
                  << " lowest_pair_power_change="
                  << (finest_pair_power - coarsest_pair_power) /
                     coarsest_pair_power << "\n";

        std::cout << "circular_vector_fem_modes=" << finest_mode_count << "\n";
        std::cout << "maximum_rotation_covariance_error=" << maximum_rotation_error << "\n";
        std::cout << "circular_localized_strike_projection: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "circular_localized_strike_projection: FAIL: " << e.what() << "\n";
        return 1;
    }
}
