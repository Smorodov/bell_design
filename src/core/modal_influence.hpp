#pragma once

#include "bell_types.h"
#include "modal_types.hpp"

namespace BellOptimizer {

// Supplies a neutral generalized influence field to the modal solver. The FEM
// solver only consumes the returned vectors and projects them through its mass
// operator; it does not know whether the field represents a strike, a rotation
// proxy, a pressure load, or another excitation model.
class ModalInfluenceProvider {
public:
    virtual ~ModalInfluenceProvider() = default;
    virtual ModalInfluence make(const BellMeshData& mesh) const = 0;
    virtual const char* name() const noexcept = 0;
};


} // namespace BellOptimizer
