#include "meridional_influence.hpp"

namespace BellOptimizer::Excitation {

ModalInfluence MeridionalInfluence::make(const BellMeshData& mesh) const {
    ModalInfluence influence;
    influence.id = name();
    influence.reference_z = reference_z_;
    influence.radial.assign(mesh.nodes.size(), 0.0);
    influence.axial.assign(mesh.nodes.size(), 0.0);

    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const auto& node = mesh.nodes[i];
        influence.radial[i] = -(node.z - reference_z_);
        influence.axial[i] = node.r;
    }
    return influence;
}

} // namespace BellOptimizer::Excitation
