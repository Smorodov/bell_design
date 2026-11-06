#include "acoustics/finite_distance_radiation.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;

int main() {
    try {
        BellGeometry geometry;
        geometry.H = 658.3175709034034;
        geometry.R_mouth = 355.9758748202059;
        geometry.t_top = 22.85675347179014;
        geometry.t_bottom = 46.62524667849419;
        geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
        geometry.ellipticity = 0.9999;
        geometry.scale(1000.0 / geometry.H);

        BellMaterial material;
        BellMeshConfig mesh;
        mesh.nNodes = 17;
        const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh);

        FEM::CylindricalModalEvaluationConfig config;
        config.harmonics = {0, 1};
        config.angular_samples = 48;
        config.requested_modes = 8;
        config.max_iterations = 3000;
        config.tolerance = 1.0e-6;
        config.minimum_frequency_hz = 1.0;
        ModalData modes;
        Excitation::MeridionalInfluence influence;
        if (!FEM::evaluate_cylindrical_bell_modes(
                model, influence.make(model.mesh_si), config, modes))
            throw std::runtime_error("real FEM modal solve failed");
        if (modes.modes.empty()) throw std::runtime_error("real FEM returned no modes");

        const auto& profile = model.mesh_si.outer_profile;
        if (profile.size() != static_cast<std::size_t>(mesh.nNodes))
            throw std::runtime_error("unexpected production outer-profile size");

        // Exercise several actual FEM modes at practical listening distances.
        // The angular quadrature is refined independently for every receiver.
        // This is a convergence contract, not an assertion that the current
        // open-shell Green-kernel proxy is an exact exterior acoustic solution.
        const double distances[] = {1.0, 2.0, 3.0, 5.0};
        const double polar = 0.82;
        const double azimuth = 0.37;
        const double earSpacing = 0.18;
        double worstQuadratureError = 0.0;
        double worstMeridionalError = 0.0;
        std::size_t evaluatedModes = 0;
        bool sawNonAxisymmetricMode = false;

        for (const auto& mode : modes.modes) {
            if (mode.outer_surface_harmonics.empty()) continue;
            const bool hasNonAxisymmetric = std::any_of(
                mode.outer_surface_harmonics.begin(), mode.outer_surface_harmonics.end(),
                [](int h) { return h != 0; });
            sawNonAxisymmetricMode = sawNonAxisymmetricMode || hasNonAxisymmetric;

            for (double distance : distances) {
                const auto listener = Acoustics::makeStereoListenerGeometry(
                    distance, polar, azimuth, earSpacing, 0.0);
                const auto leftCoarse = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.left, 360);
                const auto leftFine = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.left, 720);
                const auto leftRefined = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.left, 1440);
                const auto rightCoarse = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.right, 360);
                const auto rightFine = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.right, 720);
                const auto rightRefined = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.right, 1440);
                const auto leftMeridional8 = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.left, 720, 343.0, 8);
                const auto leftMeridional16 = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.left, 720, 343.0, 16);
                const auto rightMeridional8 = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.right, 720, 343.0, 8);
                const auto rightMeridional16 = Acoustics::evaluateFiniteDistanceRadiation(
                    mode, profile, listener.right, 720, 343.0, 16);

                const auto relativeError = [](const std::complex<double>& a,
                                              const std::complex<double>& b) {
                    return std::abs(a-b) / std::max(1.0e-24, std::abs(b));
                };
                const double leftError = relativeError(leftFine, leftRefined);
                const double rightError = relativeError(rightFine, rightRefined);
                const double leftMeridionalError = relativeError(leftMeridional8, leftMeridional16);
                const double rightMeridionalError = relativeError(rightMeridional8, rightMeridional16);
                worstQuadratureError = std::max(
                    worstQuadratureError, std::max(leftError, rightError));
                worstMeridionalError = std::max(
                    worstMeridionalError, std::max(leftMeridionalError, rightMeridionalError));
                if (!(std::isfinite(leftRefined.real()) && std::isfinite(leftRefined.imag()) &&
                      std::isfinite(rightRefined.real()) && std::isfinite(rightRefined.imag()) &&
                      std::abs(leftRefined) > 1.0e-24 && std::abs(rightRefined) > 1.0e-24))
                    throw std::runtime_error("real-FEM receiver field is invalid or numerically zero");
                if (!(leftMeridionalError < 5.0e-3 && rightMeridionalError < 5.0e-3)) {
                    std::cerr << "meridional quadrature convergence failed: f=" << mode.frequency_hz
                              << " Hz, R=" << distance << " m, L=" << leftMeridionalError
                              << ", R=" << rightMeridionalError << '\n';
                    return 4;
                }
                if (!(leftError < 5.0e-3 && rightError < 5.0e-3)) {
                    std::cerr << "quadrature convergence failed: f=" << mode.frequency_hz
                              << " Hz, R=" << distance << " m, L=" << leftError
                              << ", R=" << rightError << '\n';
                    return 2;
                }
                // A 360 -> 720 refinement should not diverge catastrophically;
                // this catches aliasing while allowing the final pair to define
                // the quantitative convergence threshold above.
                if (!(relativeError(leftCoarse, leftFine) < 2.0e-2 &&
                      relativeError(rightCoarse, rightFine) < 2.0e-2))
                    throw std::runtime_error("coarse azimuth quadrature is unstable");
            }
            ++evaluatedModes;
        }
        if (evaluatedModes < 2 || !sawNonAxisymmetricMode)
            throw std::runtime_error("real FEM test did not exercise multiple/non-axisymmetric modes");

        // Retain an explicit symmetry check for an axisymmetric mode if one
        // was produced: identical geometric conditions must give L == R.
        bool checkedAxisymmetric = false;
        for (const auto& mode : modes.modes) {
            const bool axisymmetric = mode.outer_surface_harmonics.size() == 1 &&
                                      mode.outer_surface_harmonics.front() == 0;
            if (!axisymmetric) continue;
            const auto listener = Acoustics::makeStereoListenerGeometry(
                3.0, polar, azimuth, earSpacing, 0.0);
            const auto left = Acoustics::evaluateFiniteDistanceRadiation(
                mode, profile, listener.left, 720);
            const auto right = Acoustics::evaluateFiniteDistanceRadiation(
                mode, profile, listener.right, 720);
            const double symmetryError = std::abs(left-right) /
                std::max({1.0e-24, std::abs(left), std::abs(right)});
            if (!(std::isfinite(symmetryError) && symmetryError < 2.0e-8)) {
                std::cerr << "axisymmetric real-FEM L/R symmetry error="
                          << symmetryError << '\n';
                return 3;
            }
            checkedAxisymmetric = true;
            break;
        }
        if (!checkedAxisymmetric)
            std::cout << "axisymmetric mode not returned; symmetry subcheck skipped\n";

        std::cout << "finite_distance_real_fem_contract: PASS; modes="
                  << evaluatedModes << "; distances=1,2,3,5 m; max 720->1440 error="
                  << worstQuadratureError << "; max meridional 8->16 error="
                  << worstMeridionalError << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "finite_distance_real_fem_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
