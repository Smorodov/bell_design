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
    // Structural loss factor (approximately 2*zeta for light damping).
    // A value is provided only when supported by bell-bronze measurements.
    double structural_loss_factor;
    bool has_structural_loss_factor;
    const char* damping_source;
};

// Material properties are the mechanical inputs used by BellOptimizer.
// Chemical composition is intentionally outside the solver configuration.
const BellMaterialPreset* find_material_preset(const std::string& id);
const BellMaterialPreset& default_material_preset();
std::vector<std::string> material_preset_ids();

} // namespace BellOptimizer

#endif // BELL_MATERIAL_H
