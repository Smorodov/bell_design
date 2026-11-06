#pragma once

#include "directional_radiation.hpp"

#include <string>
#include <vector>

namespace BellOptimizer::Acoustics {

// Caller-supplied stable hashes bind a cache to the exact bell geometry,
// material preset and FEM/modal solver configuration.
struct AcousticFieldCacheIdentity {
    std::string geometry_fingerprint;
    std::string material_fingerprint;
    std::string solver_fingerprint;
};

struct CachedModalRadiationField {
    double frequency_hz = 0.0;
    AcousticPolarDirectionalField directional_field;
};

// Geometry/material-dependent radiation data only. No excitation force,
// clapper trajectory, contact history, modal release state, T60 or PCM is stored.
// The canonical 0..pi polar grid makes this reusable across listener positions.
struct AcousticFieldCache {
    AcousticFieldCacheIdentity identity;
    std::vector<CachedModalRadiationField> modes;
    bool valid() const;
};

AcousticFieldCache buildAcousticFieldCache(
    const std::vector<const ModalMode*>& modes,
    const std::vector<Point2D>& outerProfile,
    const AcousticFieldCacheIdentity& identity,
    unsigned polarSamples = 91);

void saveAcousticFieldCache(const std::string& path,
                            const AcousticFieldCache& cache);
AcousticFieldCache loadAcousticFieldCache(
    const std::string& path, const AcousticFieldCacheIdentity& expectedIdentity);

// Construct a lightweight receiver-query state from a loaded cache. This does
// no FEM solve and no radiation surface integration. Distance, centre angle,
// azimuth, ear spacing and yaw envelope are per-listening-session parameters.
PreparedStereoRadiation prepareStereoRadiationFromCache(
    const AcousticFieldCache& cache, double centreDistanceM,
    double centrePolarRad, double centreAzimuthRad, double earSpacingM,
    double maximumAbsYawRad = 1.5707963267948966);

} // namespace BellOptimizer::Acoustics
