#include "models/bell_model_profile.hpp"
#include <cassert>
#include <iostream>

int main() {
    BellMaterial material;

    for (const char* family : {"neutral", "european", "chinese", "orthodox"}) {
        auto provider = BellOptimizer::Models::make_seed_provider(family);
        assert(provider);
        const auto seeds = provider->make_seeds(material, 7);
        assert(!seeds.empty());
        for (const auto& g : seeds) assert(g.isValid());
    }

    assert(BellOptimizer::Models::is_supported_tuning_model("western_five_partial"));
    assert(!BellOptimizer::Models::is_supported_tuning_model("dual_tone_chinese"));

    BellGeometry g;
    g.r_knots = {0.0, 0.2, 0.5, 0.8, 1.0};
    const auto zh = BellOptimizer::Models::assess_tradition(g, "chinese_zhong");
    assert(!zh.compatible_with_axisymmetric_model);

    std::cout << "PASS: concrete bell families remain in the model layer.\n";
    return 0;
}
