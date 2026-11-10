#include "models/bell_model_profile.hpp"
#include <cmath>
#include <iostream>

int main() {
    BellGeometry g;
    g.H = 1000.0;
    g.R_mouth = 520.0;
    g.t_top = 10.0;
    g.t_bottom = 35.0;
    g.r_knots = {0.0, 0.13, 0.36, 0.69, 1.0};

    const auto eu = BellOptimizer::Models::assess_tradition(g, "european");
    const auto ru = BellOptimizer::Models::assess_tradition(g, "russian_orthodox");
    const auto zh = BellOptimizer::Models::assess_tradition(g, "chinese_zhong");

    if (!(eu.geometry_score > 40.0 && ru.geometry_score > 40.0)) return 1;
    if (zh.compatible_with_axisymmetric_model) return 2;
    if (!BellOptimizer::Models::is_supported_tuning_model("russian_orthodox_sysoi_2000_experimental")) return 4;
    if (BellOptimizer::Models::is_supported_tuning_model("unsupported_orthodox_guess")) return 5;
    if (!std::isfinite(eu.geometry_score) || !std::isfinite(ru.geometry_score)) return 3;
    std::cout << "PASS: tradition assessment and experimental Russian Orthodox tuning-model registration.\n";
    return 0;
}
