#pragma once

#include "bell_types.h"

#include <string>
#include <vector>

namespace BellOptimizer::Validation {

struct ReferenceNode3D {
    int id = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct ReferenceHex8 {
    int n[8]{};
};

struct NonAxisymmetricReferenceMesh {
    std::vector<ReferenceNode3D> nodes;
    std::vector<ReferenceHex8> elements;
    double physical_volume = 0.0;
    double computational_volume = 0.0;
};

// Validation-only 3-D reference geometry. Production FEM remains axisymmetric.
// The meridional Q4 mesh is revolved into a periodic C3D8 mesh and then mapped
// with x=xc, y=q*yc. This is intentionally kept outside the production solver.
NonAxisymmetricReferenceMesh build_reference_mesh(
    const BellMeshData& axisymmetric_mesh,
    double q,
    int azimuthal_divisions);

bool write_calculix_c3d8_frequency_input(
    const NonAxisymmetricReferenceMesh& mesh,
    const BellMaterial& material,
    const std::string& filename,
    int requested_modes,
    std::string* error = nullptr);

} // namespace BellOptimizer::Validation
