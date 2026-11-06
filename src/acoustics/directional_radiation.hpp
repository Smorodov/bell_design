#pragma once

#include "acoustic_directional_field.hpp"
#include "../core/bell_types.h"
#include "../core/modal_types.hpp"
#include "listener_geometry.hpp"

#include <complex>
#include <string>
#include <vector>

namespace BellOptimizer::Acoustics {

struct AcousticFieldCache;

// Direct azimuthal quadrature reference implementation.
std::complex<double> evaluateDirectionalRadiationDirect(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double listenerPolarRad, double listenerAzimuthRad,
    unsigned azimuthSamples = 96);

// Fourier/Bessel representation of the same far-field surface integral.
AcousticDirectionalField buildDirectionalRadiationField(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double listenerPolarRad, double observationRadiusM);

// Precomputed polar-angle family of azimuthal fields. The FEM/modal surface
// data is integrated once per polar node; arbitrary receiver directions are
// then evaluated by complex-coefficient interpolation plus azimuthal phase
// summation. This avoids repeating FEM/contact evolution and avoids rebuilding
// a surface integral for each ear or head-yaw query.
struct AcousticPolarDirectionalField {
    std::vector<double> polar_angles_rad;
    std::vector<AcousticDirectionalField> fields;

    bool valid() const;
    std::complex<double> evaluate(double polarRad, double azimuthRad) const;
    std::complex<double> evaluateRotated(double polarRad, double azimuthRad,
                                         double frameRotationRad) const;
};

AcousticPolarDirectionalField buildDirectionalRadiationPolarField(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double observationRadiusM, unsigned polarSamples = 91);

// Build only the polar interval needed by a receiver configuration. The
// endpoints are included exactly; this is the production stereo path so the
// radiation surface integral is not repeated over the entire sphere.
AcousticPolarDirectionalField buildDirectionalRadiationPolarFieldRange(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double observationRadiusM, double minimumPolarRad,
    double maximumPolarRad, unsigned polarSamples);


// Immutable, reusable far-field stereo receiver state. FEM modes and their
// modal/contact release state are prepared by the caller; this object builds
// the per-mode polar radiation families once and serves any receiver yaw in
// the declared envelope without rebuilding a surface integral.
struct StereoRadiationTransfer {
    std::complex<double> left;
    std::complex<double> right;
    double left_range_scale = 1.0;
    double right_range_scale = 1.0;
    double left_delay_s = 0.0;
    double right_delay_s = 0.0;
};

class PreparedStereoRadiation {
public:
    std::vector<StereoRadiationTransfer> query(double headYawRad,
                                               double centreAzimuthOffsetRad = 0.0) const;
    std::size_t modeCount() const { return fields_.size(); }
    double modeFrequencyHz(std::size_t index) const { return mode_frequencies_hz_.at(index); }
    std::size_t surfaceIntegralBuildCount() const { return surface_integral_build_count_; }
    std::size_t polarNodeCount() const { return surfaceIntegralBuildCount(); }

private:
    friend PreparedStereoRadiation prepareStereoRadiationFromCache(const AcousticFieldCache&, double, double, double, double, double);
    friend PreparedStereoRadiation prepareStereoRadiation(
        const std::vector<const ModalMode*>&, const std::vector<Point2D>&,
        double, double, double, double, double, double);
    std::vector<AcousticPolarDirectionalField> fields_;
    std::vector<double> mode_frequencies_hz_;
    std::size_t surface_integral_build_count_ = 0;
    double centreDistanceM_ = 0.0;
    double centrePolarRad_ = 0.0;
    double centreAzimuthRad_ = 0.0;
    double earSpacingM_ = 0.0;
    double maximumAbsYawRad_ = 0.0;
};

PreparedStereoRadiation prepareStereoRadiation(
    const std::vector<const ModalMode*>& modes,
    const std::vector<Point2D>& profile,
    double centreDistanceM, double centrePolarRad,
    double centreAzimuthRad, double earSpacingM,
    double maximumAbsYawRad = 1.5707963267948966,
    double maximumPolarStepRad = 0.03490658503988659);

} // namespace BellOptimizer::Acoustics
