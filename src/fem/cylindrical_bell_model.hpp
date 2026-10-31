#pragma once

#include "core/bell_types.h"
#include "cylindrical_global_assembly.hpp"

namespace BellOptimizer::FEM {

// Production adapter from BellGeometry/BellMaterial to the cylindrical vector FEM.
// BellGenerator remains the single owner of meridional geometry/mesh construction;
// this adapter only converts the existing production mesh/material units to SI and
// carries BellGeometry::ellipticity into the FEM assembly.
struct CylindricalBellModel {
    BellMeshData mesh_si;
    double young_modulus_pa = 0.0;
    double poisson_ratio = 0.0;
    double density_kg_m3 = 0.0;
    double ellipticity = 1.0;
};

CylindricalBellModel make_cylindrical_bell_model(
    const BellGeometry& geometry, const BellMaterial& material,
    const BellMeshConfig& mesh_config);

CylindricalGlobalSystem assemble_cylindrical_bell_system(
    const CylindricalBellModel& model,
    const std::vector<int>& harmonics = {0, 1, 2},
    unsigned angular_samples = 64, double zero_tolerance = 0.0);

} // namespace BellOptimizer::FEM
