#include "excitation/modal_impulse_response.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using BellOptimizer::Excitation::ModalImpulseTerm;
using BellOptimizer::Excitation::synthesizeModalImpulseResponse;

int main() {
    constexpr double fs = 48000.0;
    constexpr double duration = 0.25;
    std::vector<double> response;
    const std::vector<ModalImpulseTerm> one{{440.0, 2.0, 0.5, 2.0}};
    assert(synthesizeModalImpulseResponse(one, fs, duration, response));
    assert(response.size() == static_cast<std::size_t>(std::ceil(fs * duration)));
    assert(std::abs(response.front()) < 1e-15);

    // The initial displacement derivative equals J/M. Estimate it from the
    // first sample and compare with the exact impulse initial velocity.
    const double derivative = response[1] * fs;
    assert(std::abs(derivative - 0.25) < 2.0e-3);

    const std::vector<ModalImpulseTerm> negative{{440.0, 2.0, -0.5, 2.0}};
    std::vector<double> inverted;
    assert(synthesizeModalImpulseResponse(negative, fs, duration, inverted));
    for (std::size_t i = 0; i < response.size(); ++i)
        assert(std::abs(response[i] + inverted[i]) < 1e-14);

    const std::vector<ModalImpulseTerm> two{{440.0, 2.0, 0.5, 2.0},
                                             {660.0, 1.0, -0.2, 1.5}};
    std::vector<double> combined, a, b;
    assert(synthesizeModalImpulseResponse(two, fs, duration, combined));
    assert(synthesizeModalImpulseResponse({two[0]}, fs, duration, a));
    assert(synthesizeModalImpulseResponse({two[1]}, fs, duration, b));
    for (std::size_t i = 0; i < combined.size(); ++i)
        assert(std::abs(combined[i] - a[i] - b[i]) < 1e-14);

    assert(!synthesizeModalImpulseResponse({{440.0, 0.0, 1.0, 2.0}}, fs, duration, response));
    assert(!synthesizeModalImpulseResponse({{1.0, 1.0, 1.0, 0.01}}, fs, duration, response));
    assert(!synthesizeModalImpulseResponse(one, 0.0, duration, response));
    std::cout << "modal_impulse_response_contract: PASS\n";
}
