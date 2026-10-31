#include "modal_tracking.hpp"

#include <cmath>
#include <iostream>

using namespace BellOptimizer;

static ModalMode makeMode(int number, double frequency,
                          std::initializer_list<double> radial,
                          std::initializer_list<double> axial) {
    ModalMode m;
    m.mode_number = number;
    m.frequency_hz = frequency;
    m.outer_surface_radial.assign(radial.begin(), radial.end());
    m.outer_surface_axial.assign(axial.begin(), axial.end());
    return m;
}

int main() {
    ModalData reference;
    reference.modes.push_back(makeMode(1, 100.0, {1, 0, -1}, {0, 1, 0}));
    reference.modes.push_back(makeMode(2, 150.0, {0, 1, 0}, {1, 0, -1}));
    reference.modes.push_back(makeMode(3, 210.0, {1, -1, 1}, {0, 0, 0}));

    // Deliberately swap the first two frequency ranks. Eigenvector signs are
    // also reversed; identity must follow shape, not sign or frequency order.
    ModalData candidate;
    candidate.modes.push_back(makeMode(2, 145.0, {0, -1, 0}, {-1, 0, 1}));
    candidate.modes.push_back(makeMode(1, 155.0, {-1, 0, 1}, {0, -1, 0}));
    candidate.modes.push_back(makeMode(3, 205.0, {-1, 1, -1}, {0, 0, 0}));

    const auto report = trackModalFamilies(reference, candidate);
    if (!report.used_shape_tracking || report.matched_modes != 3) return 1;
    if (candidate.modes[0].mode_number != 1 ||
        candidate.modes[1].mode_number != 2 ||
        candidate.modes[2].mode_number != 3) return 2;
    if (report.minimum_similarity < 0.999 || report.mean_similarity < 0.999) return 3;
    if (report.ambiguous) return 4;

    std::cout << "modal_family_tracking: PASS"
              << " min_similarity=" << report.minimum_similarity
              << " mean_similarity=" << report.mean_similarity << "\n";
    return 0;
}
