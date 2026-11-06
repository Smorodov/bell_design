#include "acoustics/listener_geometry.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace {
constexpr double pi = 3.14159265358979323846;
bool near(double a, double b, double tolerance = 1.0e-12) {
    return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}
}

int main() {
    using BellOptimizer::Acoustics::makeStereoListenerGeometry;
    const double R = 3.0, d = 0.18, polar = 30.0 * pi / 180.0;

    const auto front = makeStereoListenerGeometry(R, polar, 0.0, d, 0.0);
    if (!near(front.left.distance_m, front.right.distance_m)) return 1;
    if (!near(front.left.polar_rad, front.right.polar_rad) ||
        !(front.left.polar_rad > polar)) return 2;
    if (!(front.left.azimuth_rad < 0.0 && front.right.azimuth_rad > 0.0)) return 3;
    if (!near(front.left.propagation_delay_s, front.right.propagation_delay_s)) return 4;
    if (!near(front.ear_angular_separation_rad,
              2.0 * std::asin(d / (2.0 * std::sqrt(R * R + 0.25 * d * d))), 1.0e-10)) return 5;

    const auto yawed = makeStereoListenerGeometry(R, polar, 0.0, d, pi / 4.0);
    if (!(yawed.left.distance_m < R && yawed.right.distance_m > R)) return 6;
    if (!(yawed.left.propagation_delay_s < 0.0 && yawed.right.propagation_delay_s > 0.0)) return 7;
    if (!(yawed.left.polar_rad < polar && yawed.right.polar_rad > polar)) return 8;

    const auto rotated = makeStereoListenerGeometry(R, polar, 0.7, d, 0.0);
    if (!near(rotated.left.distance_m, front.left.distance_m) ||
        !near(rotated.right.distance_m, front.right.distance_m)) return 9;
    if (!near(rotated.ear_angular_separation_rad, front.ear_angular_separation_rad)) return 10;

    const auto farther = makeStereoListenerGeometry(6.0, polar, 0.0, d, 0.0);
    if (!(farther.ear_angular_separation_rad < front.ear_angular_separation_rad)) return 11;

    bool rejected = false;
    try { (void)makeStereoListenerGeometry(3.0, polar, 0.0, 0.18, 2.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) return 12;

    std::cout << "stereo_listener_geometry_contract: PASS\n";
    return 0;
}
