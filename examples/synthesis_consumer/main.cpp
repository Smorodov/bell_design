#include <acoustics/bell_modal_profile.hpp>
#include <excitation/clapper_dynamics.hpp>

#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: bell_synthesis_consumer <profile-base> <geometry-fingerprint> <material-fingerprint>\n";
        return 2;
    }
    try {
        // The solver fingerprint is supplied by the host application as part
        // of its profile-selection policy; no FEM or optimizer is linked here.
        BellOptimizer::Acoustics::AcousticFieldCacheIdentity identity{
            argv[2], argv[3], "cylindrical-modal-profile-v1"};
        auto bundle = BellOptimizer::Acoustics::loadBellModalProfileBundle(argv[1], identity);
        auto modes = BellOptimizer::Acoustics::reconstructBoundaryModes(bundle.structural);
        std::cout << "Loaded " << modes.size() << " modes; profile is ready for excitation.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Bell synthesis profile error: " << e.what() << '\n';
        return 1;
    }
}
