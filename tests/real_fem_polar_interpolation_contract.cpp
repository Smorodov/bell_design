#include "acoustics/directional_radiation.hpp"
#include "acoustics/acoustic_field_cache.hpp"
#include "acoustics/bell_modal_profile.hpp"
#include <cstdio>
#include "acoustics/listener_geometry.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "excitation/meridional_influence.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::Acoustics;

namespace {
constexpr double kPi = 3.14159265358979323846;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct ErrorStats {
    double maxNormalizedComplexError = 0.0;
    double maxPhaseErrorDeg = 0.0;
    unsigned comparisons = 0;
};
}

int main() {
    try {
        ErrorStats stats;
        unsigned modesChecked = 0;
        unsigned directChecks = 0;
        unsigned geometryCases = 0;
        double minimumFrequencyHz = 1.0e300;
        double maximumFrequencyHz = 0.0;
        const double distances[] = {1.0, 3.0, 5.0};
        const double yawAngles[] = {-kPi/2.0, -kPi/4.0, 0.0, kPi/4.0, kPi/2.0};
        for (unsigned caseIndex = 0; caseIndex < 3; ++caseIndex) {
            BellGeometry geometry;
            geometry.H = 658.3175709034034;
            geometry.R_mouth = 355.9758748202059;
            geometry.t_top = 22.85675347179014;
            geometry.t_bottom = 46.62524667849419;
            geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
            geometry.ellipticity = 0.9999;
            unsigned nodeCount = 17;
            if (caseIndex == 1) nodeCount = 25; // mesh refinement on the same profile
            if (caseIndex == 2) {
                geometry.H = 700.0;
                geometry.R_mouth = 330.0;
                geometry.t_top = 24.0;
                geometry.t_bottom = 42.0;
                geometry.r_knots = {0.0, 0.11, 0.24, 0.51, 1.0};
            }
            geometry.scale(1000.0 / geometry.H);
            BellMaterial material;
            BellMeshConfig mesh;
            mesh.nNodes = nodeCount;
            const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh);
            FEM::CylindricalModalEvaluationConfig config;
            config.harmonics = {0, 1, 2};
            config.angular_samples = 48;
            config.requested_modes = 12;
            config.max_iterations = 3000;
            config.tolerance = 1.0e-6;
            config.minimum_frequency_hz = 1.0;
            ModalData modes;
            Excitation::MeridionalInfluence influence;
            require(FEM::evaluate_cylindrical_bell_modes(
                        model, influence.make(model.mesh_si), config, modes),
                    "real FEM modal solve failed");
            const auto& profile = model.mesh_si.outer_profile;
            require(!modes.modes.empty(), "real FEM returned no modes");

            std::vector<const ModalMode*> cacheModes;
            for (const auto& mode : modes.modes)
                if (!mode.outer_surface_harmonics.empty()) cacheModes.push_back(&mode);
            if (caseIndex == 0) {
                std::vector<const ModalMode*> cacheAuditModes(cacheModes.begin(),
                    cacheModes.begin() + std::min<std::size_t>(4, cacheModes.size()));
                const AcousticFieldCacheIdentity cacheIdentity{
                    "real-fem-profile-case-0", "default-material-preset", "cylindrical-modal-config-v1"};
                const auto reusableCache = buildAcousticFieldCache(
                    cacheAuditModes, profile, cacheIdentity, 91);
                require(reusableCache.valid(), "real-FEM full-polar cache failed validation");
                double cachePeak = 0.0, cacheError = 0.0;
                for (std::size_t m = 0; m < reusableCache.modes.size(); ++m) {
                    const auto& field = reusableCache.modes[m].directional_field;
                    for (double polar : {0.13, 0.51, 0.94, 1.37, 1.91, 2.43, 2.96}) {
                        const auto exact = buildDirectionalRadiationField(
                            *cacheAuditModes[m], profile, polar, 1.0);
                        for (double azimuth : {-2.1, -0.4, 0.0, 0.83, 2.2}) {
                            const auto reference = exact.evaluate(azimuth);
                            const auto value = field.evaluate(polar, azimuth);
                            cachePeak = std::max(cachePeak, std::abs(reference));
                            cacheError = std::max(cacheError, std::abs(value - reference));
                        }
                    }
                }
                const double cacheNormalizedError = cacheError / std::max(cachePeak, 1.0e-30);
                std::cout << "real_fem_full_polar_cache_modes=" << cacheAuditModes.size()
                          << " interpolation_error_over_peak=" << cacheNormalizedError << '\n';
                require(cacheNormalizedError < 0.01,
                        "real-FEM canonical full-polar cache exceeds 1% interpolation error");

                ModalData structuralSubset;
                for (std::size_t i = 0; i < std::min<std::size_t>(4, cacheAuditModes.size()); ++i)
                    structuralSubset.modes.push_back(*cacheAuditModes[i]);
                const AcousticFieldCacheIdentity profileIdentity{
                    "real-fem-modal-profile-case-0", "default-material-preset", "cylindrical-modal-config-v1"};
                const auto modalProfile = buildBellModalProfile(
                    structuralSubset, model.mesh_si.inner_profile, profile, profileIdentity);
                const std::string bundleBase = "real_fem_modal_profile_contract";
                saveBellModalProfileBundle(bundleBase, modalProfile, 91);
                const auto loadedBundle = loadBellModalProfileBundle(bundleBase, profileIdentity);
                const auto restoredModes = reconstructBoundaryModes(loadedBundle.structural);
                require(restoredModes.size() == structuralSubset.modes.size() &&
                        loadedBundle.radiation.modes.size() == structuralSubset.modes.size(),
                        "real-FEM modal profile bundle changed the mode count");
                for (std::size_t i = 0; i < restoredModes.size(); ++i) {
                    require(restoredModes[i].inner_surface_fourier_radial ==
                                structuralSubset.modes[i].inner_surface_fourier_radial &&
                            restoredModes[i].outer_surface_fourier_axial ==
                                structuralSubset.modes[i].outer_surface_fourier_axial,
                            "real-FEM modal profile bundle changed a boundary trace");
                }
                std::remove((bundleBase + ".bmodal").c_str());
                std::remove((bundleBase + ".acfield").c_str());
            }
            const auto preparedStereo = prepareStereoRadiation(
                cacheModes, profile, 3.0, 0.82, 0.37, 0.18);
            require(preparedStereo.modeCount() == cacheModes.size(),
                    "prepared stereo cache lost one or more real FEM modes");
            const std::size_t preparedSurfaceBuilds = preparedStereo.surfaceIntegralBuildCount();
            require(preparedSurfaceBuilds >= preparedStereo.modeCount(),
                    "prepared stereo cache reports no polar surface-integral work");
            for (double azimuthOffset : {0.0, 0.37, -1.1}) {
              for (double yaw : yawAngles) {
                const auto transfers = preparedStereo.query(yaw, azimuthOffset);
                require(transfers.size() == cacheModes.size(),
                        "listener query changed the prepared mode count");
                const auto ears = makeStereoListenerGeometry(3.0, 0.82, 0.37 + azimuthOffset, 0.18, yaw);
                for (std::size_t i = 0; i < transfers.size(); ++i) {
                    const auto exactLeft = buildDirectionalRadiationField(
                        *cacheModes[i], profile, ears.left.polar_rad, 3.0).evaluate(ears.left.azimuth_rad);
                    const auto exactRight = buildDirectionalRadiationField(
                        *cacheModes[i], profile, ears.right.polar_rad, 3.0).evaluate(ears.right.azimuth_rad);
                    const double peak = std::max({1.0e-30, std::abs(exactLeft), std::abs(exactRight)});
                    require(std::abs(transfers[i].left - exactLeft) / peak < 0.01 &&
                            std::abs(transfers[i].right - exactRight) / peak < 0.01,
                            "prepared stereo listener query exceeds 1% exact-polar error");
                    require(std::isfinite(transfers[i].left_range_scale) &&
                            std::isfinite(transfers[i].right_range_scale) &&
                            std::isfinite(transfers[i].left_delay_s) &&
                            std::isfinite(transfers[i].right_delay_s),
                            "prepared stereo listener query returned invalid receiver geometry");
                }
                require(preparedStereo.surfaceIntegralBuildCount() == preparedSurfaceBuilds,
                        "listener yaw/azimuth query rebuilt a polar surface integral");
              }
            }
            bool rejectedOutsideYaw = false;
            try { (void)preparedStereo.query(0.5 * kPi + 1.0e-4); }
            catch (const std::out_of_range&) { rejectedOutsideYaw = true; }
            require(rejectedOutsideYaw, "prepared stereo cache accepted yaw outside its declared envelope");

            unsigned caseModes = 0;
            for (const auto& mode : modes.modes) {
                if (mode.outer_surface_harmonics.empty()) continue;
                minimumFrequencyHz = std::min(minimumFrequencyHz, mode.frequency_hz);
                maximumFrequencyHz = std::max(maximumFrequencyHz, mode.frequency_hz);
                double modePeak = 0.0;
                double modeMaxComplexError = 0.0;
                std::vector<std::pair<std::complex<double>, std::complex<double>>> modeSamples;
                for (double distance : distances) {
                    const auto yawA = makeStereoListenerGeometry(distance, 0.82, 0.37, 0.18, -kPi/2.0);
                    const auto yawB = makeStereoListenerGeometry(distance, 0.82, 0.37, 0.18,  kPi/2.0);
                    const double pMin = std::min({yawA.left.polar_rad, yawA.right.polar_rad,
                                                  yawB.left.polar_rad, yawB.right.polar_rad});
                    const double pMax = std::max({yawA.left.polar_rad, yawA.right.polar_rad,
                                                  yawB.left.polar_rad, yawB.right.polar_rad});
                    const unsigned nodes = std::max(2U, static_cast<unsigned>(
                        std::ceil((pMax - pMin) / (2.0 * kPi / 180.0))) + 1U);
                    const auto field = buildDirectionalRadiationPolarFieldRange(
                        mode, profile, distance, pMin, pMax, nodes);
                    for (double yaw : yawAngles) {
                        const auto ears = makeStereoListenerGeometry(distance, 0.82, 0.37, 0.18, yaw);
                        for (const ReceiverPoint* ear : {&ears.left, &ears.right}) {
                            const auto interpolated = field.evaluate(ear->polar_rad, ear->azimuth_rad);
                            const auto exactField = buildDirectionalRadiationField(
                                mode, profile, ear->polar_rad, distance);
                            const auto exact = exactField.evaluate(ear->azimuth_rad);
                            modePeak = std::max(modePeak, std::abs(exact));
                            modeMaxComplexError = std::max(modeMaxComplexError,
                                                           std::abs(interpolated - exact));
                            modeSamples.emplace_back(interpolated, exact);
                            ++stats.comparisons;
                        }
                    }
                    // Independent direct azimuthal quadrature spot check at
                    // the centre yaw, separating interpolation from Bessel expansion.
                    const auto middle = makeStereoListenerGeometry(distance, 0.82, 0.37, 0.18, 0.0);
                    const auto direct = evaluateDirectionalRadiationDirect(
                        mode, profile, middle.left.polar_rad, middle.left.azimuth_rad, 4096);
                    const auto exact = buildDirectionalRadiationField(
                        mode, profile, middle.left.polar_rad, distance).evaluate(middle.left.azimuth_rad);
                    const double directError = std::abs(direct - exact) /
                        std::max({1.0e-30, std::abs(direct), std::abs(exact)});
                    require(std::isfinite(directError) && directError < 2.0e-5,
                            "real FEM exact-polar Fourier/Bessel field disagrees with direct quadrature");
                    ++directChecks;
                }
                const double modeNormalizedError = modeMaxComplexError / std::max(modePeak, 1.0e-30);
                stats.maxNormalizedComplexError = std::max(stats.maxNormalizedComplexError,
                                                            modeNormalizedError);
                for (const auto& sample : modeSamples) {
                    if (std::abs(sample.second) < 0.1 * std::max(modePeak, 1.0e-30)) continue;
                    const double phase = std::abs(std::arg(sample.first * std::conj(sample.second))) * 180.0 / kPi;
                    stats.maxPhaseErrorDeg = std::max(stats.maxPhaseErrorDeg, phase);
                }
                ++caseModes;
            }
            require(caseModes >= 4, "too few real FEM modes were exercised in a geometry case");
            modesChecked += caseModes;
            ++geometryCases;
            std::cout << "real_fem_case=" << caseIndex << " mesh_nodes=" << nodeCount
                      << " modes=" << caseModes << '\n';
        }
        require(geometryCases == 3, "not all real-FEM geometry/mesh cases ran");
        const double normalizedError = stats.maxNormalizedComplexError;
        std::cout << "real_fem_geometry_cases=" << geometryCases << " real_fem_modes=" << modesChecked
                  << " frequency_range_hz=" << minimumFrequencyHz << ".." << maximumFrequencyHz
                  << " distances=1,3,5m yaw_samples=5 ears=2 comparisons=" << stats.comparisons
                  << " direct_quadrature_checks=" << directChecks
                  << " worst_per_mode_complex_error_over_peak=" << normalizedError
                  << " max_phase_error_deg_above_10pct=" << stats.maxPhaseErrorDeg << '\n';
        require(normalizedError <= 0.01,
                "real FEM polar interpolation exceeds 1% per-mode peak-normalized complex error");
        require(stats.maxPhaseErrorDeg <= 1.0,
                "real FEM polar interpolation exceeds 1 degree phase error above 10% of each mode peak");
        std::cout << "real_fem_polar_interpolation_contract: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "real_fem_polar_interpolation_contract: FAIL: " << error.what() << '\n';
        return 1;
    }
}
