#include "excitation/clapper_modal_state_provider.hpp"

#include <iostream>
#include <vector>

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::Excitation;

    ClapperModalEvaluationConfig config;
    ClapperModalEvaluationResult result;
    result.first_contact_z_m = 123.0; // The API must clear stale output on rejection.

    const std::vector<Point2D> profile{{0.10, 0.0}, {0.20, 1.0}};
    const std::vector<ModalMode> insufficient_modes(4);
    if (evaluateClapperModalStateFromCandidate(config, profile, insufficient_modes, result)) {
        std::cerr << "candidate path must reject fewer than five target modes\n";
        return 1;
    }
    if (result.first_contact_z_m != 0.0 || !result.modal_pressure_magnitude.empty()) {
        std::cerr << "candidate path must clear output on rejection\n";
        return 1;
    }

    const std::vector<ModalMode> five_modes(5);
    if (evaluateClapperModalStateFromCandidate(config, {}, five_modes, result)) {
        std::cerr << "candidate path must reject an incomplete inner profile\n";
        return 1;
    }

    // This entry point accepts only already-computed modal data and a boundary
    // profile. Its type-level contract has no geometry/material/FEM arguments,
    // so candidate integration cannot accidentally trigger another FEM solve.
    std::cout << "clapper_candidate_data_reuse_contract: PASS\n";
    return 0;
}
