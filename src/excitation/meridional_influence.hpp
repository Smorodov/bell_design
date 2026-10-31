#pragma once

#include "../core/modal_influence.hpp"

namespace BellOptimizer::Excitation {

// Current production meridional influence field. It is a mathematical
// generalized influence, not a claim of 3-D rigid rotation.
class MeridionalInfluence final : public ModalInfluenceProvider {
public:
    explicit MeridionalInfluence(double reference_z = 0.0) : reference_z_(reference_z) {}

    ModalInfluence make(const BellMeshData& mesh) const override;
    const char* name() const noexcept override { return "meridional_q_minus_z_r"; }

private:
    double reference_z_ = 0.0;
};

} // namespace BellOptimizer::Excitation
