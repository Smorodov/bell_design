#include "../src/core/bell_types.h"
#include <cmath>
#include <iostream>

int main() {
    BellGeometry circular;
    circular.r_knots = {0.0, 0.25, 0.6, 1.0};
    if (std::abs(circular.ellipticity - 1.0) > 1e-15 || !circular.isValid()) {
        std::cerr << "default geometry must be circular and valid\n";
        return 1;
    }
    EvolutionParams defaults;
    if (defaults.optimize_ellipticity) {
        std::cerr << "ellipticity search must be opt-in; circular optimization freezes q\n";
        return 2;
    }
    circular.ellipticity = 0.85;
    if (!circular.isValid()) {
        std::cerr << "valid elliptical geometry rejected\n";
        return 3;
    }
    circular.ellipticity = 0.69;
    if (circular.isValid()) {
        std::cerr << "ellipticity below supported bound accepted\n";
        return 4;
    }
    circular.ellipticity = 1.01;
    if (circular.isValid()) {
        std::cerr << "ellipticity above circular limit accepted\n";
        return 5;
    }
    std::cout << "ellipticity_search_policy_contract: PASS\n";
    return 0;
}
