#pragma once

#include "bell_types.h"
#include <string>

namespace BellOptimizer {

class CCXInputWriter {
public:
    static bool write(const BellMeshData& mesh,
                      const BellMaterial& material,
                      const std::string& filename,
                      int requested_modes = 10);
};

} // namespace BellOptimizer
