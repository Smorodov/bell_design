#include "excitation/clapper_modal_evaluator.hpp"
#include <iostream>
#include <vector>

int main() {
    using namespace BellOptimizer;
    using namespace BellOptimizer::Excitation;
    ClapperModalEvaluationConfig config;
    ClapperModalEvaluationResult result;
    std::vector<Point2D> profile{{0.1, 0.0}, {0.2, 1.0}};
    std::vector<ModalMode> modes;
    if (evaluateClapperModalState(config, profile, modes, result)) {
        std::cerr << "empty modal state must be rejected\n";
        return 1;
    }
    config.top_z_m = 1.0;
    config.mouth_z_m = 1.0;
    if (evaluateClapperModalState(config, profile, modes, result)) {
        std::cerr << "degenerate axial reference must be rejected\n";
        return 1;
    }
    config.top_z_m = 0.0;
    config.mouth_z_m = 1.0;
    config.admissible_strike_fraction_min = 0.8;
    config.admissible_strike_fraction_max = 0.2;
    if (evaluateClapperModalState(config, profile, modes, result)) {
        std::cerr << "inverted strike interval must be rejected\n";
        return 1;
    }
    std::cout << "clapper_modal_evaluator_contract: PASS\n";
    return 0;
}
