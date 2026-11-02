#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double PI = 3.1415926535897932384626433832795;
void req(bool v, const char* m) { if (!v) throw std::runtime_error(m); }

// Diagnostic only: unit-area half-sine force pulse.
// This is not a claim that a Chinese-bell striker follows this law.
double pulse(double t, double T) {
    if (t < 0.0 || t > T) return 0.0;
    return (PI / (2.0 * T)) * std::sin(PI * t / T);
}

double area(double T, int n) {
    const double dt = T / static_cast<double>(n);
    double sum = 0.5 * (pulse(0.0, T) + pulse(T, T));
    for (int i = 1; i < n; ++i) sum += pulse(i * dt, T);
    return sum * dt;
}

double spectrum(double T, double f, int n) {
    const double dt = T / static_cast<double>(n);
    double re = 0.0, im = 0.0;
    auto add = [&](double t, double w) {
        const double p = pulse(t, T);
        const double a = 2.0 * PI * f * t;
        re += w * p * std::cos(a);
        im -= w * p * std::sin(a);
    };
    add(0.0, 0.5); add(T, 0.5);
    for (int i = 1; i < n; ++i) add(i * dt, 1.0);
    return std::hypot(re, im) * dt;
}
}

int main() {
    try {
        const std::vector<double> durations = {0.0002, 0.0005, 0.0010, 0.0020};
        const double probe_hz = 2000.0;
        double previous = 1e100;
        std::cout << std::setprecision(12)
                  << "duration_s,unit_area,normalized_spectrum_at_2000Hz\n";
        for (double T : durations) {
            const double a = area(T, 4000);
            const double s = spectrum(T, probe_hz, 4000);
            std::cout << T << ',' << a << ',' << s << '\n';
            req(std::abs(a - 1.0) < 2e-7, "temporal pulse is not unit-area");
            req(std::isfinite(s) && s >= 0.0, "non-finite temporal spectrum");
            // For this probe frequency the selected durations lie on the first
            // monotone branch. The check is mathematical, not a physical gate.
            req(s < previous, "high-frequency content did not decrease with duration");
            previous = s;
        }
        std::cout << "chinese_bell_temporal_impulse_audit: PASS\n"
                  << "status=diagnostic_kernel_property_no_physical_calibration\n"
                  << "note=half-sine pulse is an audit kernel only; no production parameter added\n";
    } catch (const std::exception& e) {
        std::cerr << "chinese_bell_temporal_impulse_audit: FAIL: " << e.what() << '\n';
        return 1;
    }
}
