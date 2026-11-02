#include "bell_material.h"

#include <cmath>
#include <iostream>

int main() {
    const auto& default_preset = BellOptimizer::default_material_preset();
    if (std::string(default_preset.id) != "bell_bronze_reference") return 1;

    const auto* bronze = BellOptimizer::find_material_preset("bell_bronze_reference");
    if (!bronze) return 2;
    if (!(bronze->material.E > 0.0 &&
          bronze->material.nu > -1.0 && bronze->material.nu < 0.5 &&
          bronze->material.rho > 0.0)) return 3;

    const auto* aluminum_bronze =
        BellOptimizer::find_material_preset("bronze_aluminum_reference");
    if (!aluminum_bronze) return 4;
    if (std::abs(aluminum_bronze->material.E - 1.10e5) > 1e-12) return 5;
    if (std::abs(aluminum_bronze->material.rho - 7.40e-9) > 1e-18) return 6;

    if (BellOptimizer::find_material_preset("Copper") != nullptr) return 7;

    const auto ids = BellOptimizer::material_preset_ids();
    if (ids.size() != 2) return 8;

    std::cout << "material presets: " << ids.size() << "\n";
    return 0;
}
