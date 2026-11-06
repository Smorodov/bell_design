#include "acoustics/bell_modal_profile.hpp"
#include "excitation/clapper_dynamics.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using namespace BellOptimizer;
    bool rejectedMissingProfile = false;
    try {
        (void)Acoustics::loadBellModalProfileBundle(
            "definitely_missing_bell_profile",
            {"geometry", "material", "solver"});
    } catch (const std::exception&) {
        rejectedMissingProfile = true;
    }
    if (!rejectedMissingProfile) {
        std::cerr << "SDK profile loader did not report a missing profile\n";
        return 1;
    }

    Excitation::ClapperDynamicsConfig config;
    Excitation::ClapperDynamicsResult result;
    if (Excitation::simulateInternalClapperContact(config, {}, {}, result)) {
        std::cerr << "SDK accepted an empty structural profile\n";
        return 1;
    }
    std::cout << "Standalone synthesis SDK link contract passed\n";
    return 0;
}
