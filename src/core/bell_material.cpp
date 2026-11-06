#include "bell_material.h"

#include <stdexcept>

namespace BellOptimizer {
namespace {

constexpr BellMaterialPreset kPresets[] = {
    {
        "bell_bronze_reference",
        "Current BellOptimizer production reference bronze; not a historical alloy certification.",
        {1.10e4, 0.34, 8.60e-9},
        1.8e-4, true,
        "Literature reference: mean modal material damping zeta=0.9e-4; converted to loss factor eta~=2*zeta. Roozen-Kroon (1992), DOI 10.6100/IR373470. Not alloy-batch calibration."
    },
    {
        "bronze_aluminum_reference",
        "Aluminum bronze reference from the supplied CalculiX material file.",
        {1.10e5, 0.30, 7.40e-9},
        0.0, false,
        "No bell-specific measured damping value assigned; use fixed T60 or provide an explicit sensitivity override."
    }
};

constexpr std::size_t kPresetCount = sizeof(kPresets) / sizeof(kPresets[0]);

} // namespace

const BellMaterialPreset* find_material_preset(const std::string& id) {
    for (std::size_t i = 0; i < kPresetCount; ++i) {
        if (id == kPresets[i].id) return &kPresets[i];
    }
    return nullptr;
}

const BellMaterialPreset& default_material_preset() {
    return kPresets[0];
}

std::vector<std::string> material_preset_ids() {
    std::vector<std::string> ids;
    ids.reserve(kPresetCount);
    for (const auto& preset : kPresets) ids.emplace_back(preset.id);
    return ids;
}

} // namespace BellOptimizer
