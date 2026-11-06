#include "excitation/modal_decay.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace BellOptimizer::Excitation;
int main() {
    ModalDecayParameters low, high, fixedA, fixedB;
    constexpr double eta = 0.002;
    assert(calculateModalDecay(500.0, ModalDecayModel::StructuralLossFactor, 5.0, eta, low));
    assert(calculateModalDecay(1000.0, ModalDecayModel::StructuralLossFactor, 5.0, eta, high));
    assert(std::abs(high.alpha_per_second / low.alpha_per_second - 2.0) < 1e-12);
    assert(std::abs(low.t60_seconds / high.t60_seconds - 2.0) < 1e-12);
    assert(std::abs(std::exp(-low.alpha_per_second * low.t60_seconds) - 0.001) < 1e-14);
    assert(calculateModalDecay(500.0, ModalDecayModel::FixedT60, 4.0, eta, fixedA));
    assert(calculateModalDecay(1000.0, ModalDecayModel::FixedT60, 4.0, eta, fixedB));
    assert(std::abs(fixedA.alpha_per_second - fixedB.alpha_per_second) < 1e-15);
    assert(std::abs(fixedA.t60_seconds - 4.0) < 1e-12);
    assert(!calculateModalDecay(0.0, ModalDecayModel::StructuralLossFactor, 4.0, eta, low));
    assert(!calculateModalDecay(500.0, ModalDecayModel::StructuralLossFactor, 4.0, 0.0, low));
    assert(!calculateModalDecay(500.0, ModalDecayModel::StructuralLossFactor, 4.0, 0.1, low));
    assert(!calculateModalDecay(500.0, ModalDecayModel::FixedT60, 0.0, eta, low));
    std::cout << "modal_decay_contract: PASS\n";
}
