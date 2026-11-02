#include "modal_impulse_response.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Excitation {
namespace {
constexpr double PI = 3.1415926535897932384626433832795;
constexpr double LOG_1000 = 6.9077552789821370521;
}

bool synthesizeModalImpulseResponse(const std::vector<ModalImpulseTerm>& terms,
                                    double sample_rate,
                                    double duration_seconds,
                                    std::vector<double>& samples) {
    samples.clear();
    if (!std::isfinite(sample_rate) || !(sample_rate > 0.0) ||
        !std::isfinite(duration_seconds) || !(duration_seconds > 0.0) ||
        terms.empty()) return false;
    const double count_real = std::ceil(sample_rate * duration_seconds);
    if (!std::isfinite(count_real) || count_real < 1.0 ||
        count_real > static_cast<double>(std::numeric_limits<std::size_t>::max()))
        return false;

    struct Prepared { double impulse_over_mass_over_wd, alpha, wd; };
    std::vector<Prepared> prepared;
    prepared.reserve(terms.size());
    for (const auto& term : terms) {
        if (!std::isfinite(term.frequency_hz) || !(term.frequency_hz > 0.0) ||
            !std::isfinite(term.modal_mass) || !(term.modal_mass > 0.0) ||
            !std::isfinite(term.generalized_impulse) ||
            !std::isfinite(term.t60_seconds) || !(term.t60_seconds > 0.0))
            return false;
        const double omega = 2.0 * PI * term.frequency_hz;
        const double alpha = LOG_1000 / term.t60_seconds;
        if (!(alpha < omega)) return false; // declared model is underdamped only
        const double wd = std::sqrt((omega - alpha) * (omega + alpha));
        if (!(wd > 0.0) || !std::isfinite(wd)) return false;
        prepared.push_back({term.generalized_impulse / (term.modal_mass * wd), alpha, wd});
    }

    const std::size_t count = static_cast<std::size_t>(count_real);
    samples.assign(count, 0.0);
    for (std::size_t n = 0; n < count; ++n) {
        const double t = static_cast<double>(n) / sample_rate;
        double value = 0.0;
        for (const auto& term : prepared)
            value += term.impulse_over_mass_over_wd *
                     std::exp(-term.alpha * t) * std::sin(term.wd * t);
        if (!std::isfinite(value)) { samples.clear(); return false; }
        samples[n] = value;
    }
    return true;
}

} // namespace BellOptimizer::Excitation
