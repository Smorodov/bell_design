#pragma once

#include "acoustic_field_cache.hpp"
#include "../core/modal_types.hpp"

#include <string>
#include <vector>

namespace BellOptimizer::Acoustics {

// Reusable structural boundary state for contact-driven modal synthesis.
// This is deliberately not a full FEM restart: stiffness/mass matrices and
// interior eigenvectors are omitted because the synthesis path consumes the
// inner/outer modal boundary traces only.
struct CachedStructuralMode {
    int mode_number = 0;
    double eigenvalue = 0.0;
    double frequency_hz = 0.0;
    double modal_mass = 0.0;
    double decay_time = 0.0;
    std::vector<int> inner_harmonics;
    std::vector<double> inner_radial, inner_azimuthal, inner_axial;
    std::vector<int> outer_harmonics;
    std::vector<double> outer_radial, outer_azimuthal, outer_axial;
};

struct BellModalProfile {
    AcousticFieldCacheIdentity identity;
    std::vector<Point2D> inner_profile;
    std::vector<Point2D> outer_profile;
    std::vector<CachedStructuralMode> modes;
    bool valid() const;
};

BellModalProfile buildBellModalProfile(
    const ModalData& modalData,
    const std::vector<Point2D>& innerProfile,
    const std::vector<Point2D>& outerProfile,
    const AcousticFieldCacheIdentity& identity);

void saveBellModalProfile(const std::string& path, const BellModalProfile& profile);
BellModalProfile loadBellModalProfile(
    const std::string& path, const AcousticFieldCacheIdentity& expectedIdentity);

// Reconstruct the exact boundary-only ModalMode data needed by the clapper
// contact evaluator and directional radiation builder. Excitation state is
// intentionally absent and must be calculated for each strike.
std::vector<ModalMode> reconstructBoundaryModes(const BellModalProfile& profile);

struct BellModalProfileBundle {
    BellModalProfile structural;
    AcousticFieldCache radiation;
};

// A logical two-file bundle: <base>.bmodal stores inner/outer modal traces;
// <base>.acfield stores reusable directional radiation fields.
void saveBellModalProfileBundle(const std::string& basePath,
                                const BellModalProfile& profile,
                                unsigned polarSamples = 91);
BellModalProfileBundle loadBellModalProfileBundle(
    const std::string& basePath, const AcousticFieldCacheIdentity& expectedIdentity);

} // namespace BellOptimizer::Acoustics
