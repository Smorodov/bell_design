#pragma once

#include "core/bell_types.h"
#include <cstddef>
#include <vector>

namespace BellOptimizer::FEM {

// Non-owning description of one full cylindrical-FEM modal field. The mesh,
// harmonic list and mode vector must outlive each comparison call.
struct CylindricalModalField {
    const BellMeshData* mesh = nullptr;
    const std::vector<int>* harmonics = nullptr;
    const std::vector<double>* mode_shape = nullptr;
    double ellipticity = 1.0;
    double density_kg_m3 = 0.0;
};

struct CylindricalModalFieldOverlap {
    bool valid = false;
    double normalized_overlap = 0.0;
    double common_volume_fraction = 0.0;
    double cross_mass = 0.0;
    double mass_norm_a = 0.0;
    double mass_norm_b = 0.0;
};

// Integrates rho * u_a . u_b over the common physical meridional domain and
// the full azimuth. Q4 fields are evaluated by their own element shape
// functions; different meshes and nested Fourier bases are supported.
CylindricalModalFieldOverlap compare_cylindrical_modal_fields(
    const CylindricalModalField& a, const CylindricalModalField& b,
    unsigned angular_samples = 64, double minimum_common_volume_fraction = 0.98);

// Minimum principal overlap for one- or two-dimensional modal families.
// The optional output reports the minimum common-domain volume fraction over
// all field pairs. A low-coverage comparison returns zero rather than silently
// comparing modes on materially different domains.
double cylindrical_modal_subspace_similarity(
    const std::vector<CylindricalModalField>& a,
    const std::vector<CylindricalModalField>& b,
    unsigned angular_samples = 64,
    double minimum_common_volume_fraction = 0.98,
    double* common_volume_fraction = nullptr);

} // namespace BellOptimizer::FEM
