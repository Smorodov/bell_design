#ifndef BELL_MATERIAL_H
#define BELL_MATERIAL_H

#include "bell_types.h"

#include <string>
#include <vector>

namespace BellOptimizer {

struct BellMaterialPreset {
    const char* id;
    const char* description;
    BellMaterial material;
};

// Material properties are the mechanical inputs used by BellOptimizer.
// Chemical composition is intentionally outside the solver configuration.
const BellMaterialPreset* find_material_preset(const std::string& id);
const BellMaterialPreset& default_material_preset();
std::vector<std::string> material_preset_ids();

} // namespace BellOptimizer

#endif // BELL_MATERIAL_H
