#pragma once

#include <vector>

namespace BellOptimizer::Excitation {

// One mass-normalized or explicitly mass-valued modal impulse term.
// generalized_impulse is the time integral of generalized force (N s in SI).
struct ModalImpulseTerm {
    double frequency_hz = 0.0;
    double modal_mass = 1.0;
    double generalized_impulse = 0.0;
    double t60_seconds = 2.0;
};

// Sums linear underdamped modal displacement responses. The amplitude is
// q_k(t)=J_k/(M_k*wd_k)*exp(-alpha_k*t)*sin(wd_k*t),
// alpha=ln(1000)/T60, wd=sqrt((2*pi*f)^2-alpha^2).
// This is an impulse-response model; it does not model finite striker contact.
bool synthesizeModalImpulseResponse(const std::vector<ModalImpulseTerm>& terms,
                                    double sample_rate,
                                    double duration_seconds,
                                    std::vector<double>& samples);

} // namespace BellOptimizer::Excitation
