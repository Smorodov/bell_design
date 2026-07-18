// =============================================================================
// bell_constants.hpp - ФИЗИЧЕСКИЕ КОНСТАНТЫ И ПАРАМЕТРЫ
// =============================================================================
#pragma once

#include <cmath>

namespace BellOptimizer {

struct BellConstants {
    static constexpr double EPS = 1e-10;
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double TWO_PI = 2.0 * PI;
    
    static constexpr double MIN_FREQUENCY = 20.0;
    static constexpr double MAX_FREQUENCY = 20000.0;
    static constexpr double MIN_AMPLITUDE = 0.0001;
    static constexpr double CLIP_LEVEL = 0.95;
    static constexpr double ATTACK_DURATION = 0.01;
};

}