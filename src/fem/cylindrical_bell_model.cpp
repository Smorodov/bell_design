#include "cylindrical_bell_model.hpp"
#include "core/bell_generator.h"
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::FEM {
namespace {
BellMeshData to_si_mesh(const BellMeshData& source) {
    BellMeshData out = source;
    constexpr double mm_to_m = 1.0e-3;
    for (auto& node : out.nodes) {
        node.r *= mm_to_m;
        node.z *= mm_to_m;
    }
    for (auto& p : out.outer_profile) {
        p.r *= mm_to_m;
        p.z *= mm_to_m;
    }
    for (auto& p : out.inner_profile) {
        p.r *= mm_to_m;
        p.z *= mm_to_m;
    }
    out.total_arc_length *= mm_to_m;
    out.min_thickness *= mm_to_m;
    out.max_thickness *= mm_to_m;
    return out;
}
}

CylindricalBellModel make_cylindrical_bell_model(
    const BellGeometry& geometry, const BellMaterial& material,
    const BellMeshConfig& mesh_config) {
    if (!geometry.isValid())
        throw std::invalid_argument("cylindrical bell model: invalid BellGeometry");
    if (!std::isfinite(material.E) || material.E <= 0.0 ||
        !std::isfinite(material.nu) || material.nu <= -1.0 || material.nu >= 0.5 ||
        !std::isfinite(material.rho) || material.rho <= 0.0)
        throw std::invalid_argument("cylindrical bell model: invalid BellMaterial");

    BellGenerator generator(material, mesh_config);
    generator.setGeometry(geometry);
    const auto validation = generator.validateGeometry();
    if (!validation.passed)
        throw std::invalid_argument("cylindrical bell model: BellGeometry validation failed: " + validation.reason);
    if (!generator.buildMesh())
        throw std::runtime_error("cylindrical bell model: production BellGenerator failed to build mesh");

    CylindricalBellModel model;
    model.mesh_si = to_si_mesh(generator.getMeshData());
    // BellMaterial uses the historical production units MPa and tonne/mm^3.
    model.young_modulus_pa = material.E * 1.0e6;
    model.poisson_ratio = material.nu;
    model.density_kg_m3 = material.rho * 1.0e12;
    model.ellipticity = geometry.ellipticity;
    return model;
}

CylindricalGlobalSystem assemble_cylindrical_bell_system(
    const CylindricalBellModel& model, const std::vector<int>& harmonics,
    unsigned angular_samples, double zero_tolerance) {
    if (model.mesh_si.nodes.empty() || model.mesh_si.elements.empty())
        throw std::invalid_argument("cylindrical bell model: empty SI mesh");
    return assemble_cylindrical_global_system(
        model.mesh_si, model.young_modulus_pa, model.poisson_ratio,
        model.density_kg_m3, harmonics, angular_samples, zero_tolerance,
        model.ellipticity);
}

} // namespace BellOptimizer::FEM
