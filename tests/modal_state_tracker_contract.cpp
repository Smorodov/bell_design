#include "core/modal_state_tracker.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using namespace BellOptimizer;

static ModalMode mode(int id, double f, int samples, int shape) {
    ModalMode m;
    m.mode_number = id;
    m.frequency_hz = f;
    for (int i = 0; i < samples; ++i) {
        const double x = static_cast<double>(i) / static_cast<double>(samples - 1);
        if (shape == 0) {
            m.outer_surface_radial.push_back(std::sin(3.141592653589793 * x));
            m.outer_surface_axial.push_back(0.25 * std::sin(2.0 * 3.141592653589793 * x));
        } else {
            m.outer_surface_radial.push_back(0.2 * std::sin(2.0 * 3.141592653589793 * x));
            m.outer_surface_axial.push_back(std::sin(3.141592653589793 * x));
        }
    }
    return m;
}

static ModalMode mixModes(const ModalMode& a, const ModalMode& b, double c, double s,
                          double frequency, int id) {
    ModalMode m;
    m.mode_number = id;
    m.frequency_hz = frequency;
    for (std::size_t i = 0; i < a.outer_surface_radial.size(); ++i) {
        m.outer_surface_radial.push_back(c * a.outer_surface_radial[i] + s * b.outer_surface_radial[i]);
        m.outer_surface_axial.push_back(c * a.outer_surface_axial[i] + s * b.outer_surface_axial[i]);
    }
    return m;
}

int main() {
    ModalData initial;
    initial.modes = {mode(1, 100.0, 9, 0), mode(2, 150.0, 9, 1)};
    ModalStateTracker tracker;
    if (!tracker.initialize(initial, 0, 1, 0.0)) return 1;

    // Different boundary sample count and swapped solver frequency ordering.
    // Eigenvector signs are arbitrary; the tracker must retain physical identity.
    ModalData next;
    auto second = mode(2, 151.8, 17, 1);
    for (auto& v : second.outer_surface_radial) v = -v;
    for (auto& v : second.outer_surface_axial) v = -v;
    auto first = mode(1, 101.2, 17, 0);
    for (auto& v : first.outer_surface_radial) v = -v;
    for (auto& v : first.outer_surface_axial) v = -v;
    next.modes = {second, first, mode(3, 240.0, 17, 0)};

    const auto report = tracker.update(next, 0.1);
    if (!report.updated || report.ambiguous) return 2;
    if (report.tracks[0].candidate_index != 1 || report.tracks[1].candidate_index != 0) return 3;
    if (report.tracks[0].shape_similarity < 0.999 || report.tracks[1].shape_similarity < 0.999) return 4;
    if (std::abs(report.tracks[0].measured_frequency_hz - 101.2) > 1e-9) return 5;

    // A near-degenerate eigensolver may rotate the basis inside the same
    // two-dimensional modal subspace. Pair tracking must remain stable even
    // when neither individual candidate resembles its original vector.
    ModalStateTracker rotated_tracker;
    if (!rotated_tracker.initialize(initial, 0, 1, 0.0)) return 6;
    const double angle = 0.63;
    ModalData rotated;
    rotated.modes = {
        mixModes(mode(1, 100.0, 17, 0), mode(2, 100.0, 17, 1), std::cos(angle), std::sin(angle), 101.2, 10),
        mixModes(mode(1, 100.0, 17, 0), mode(2, 100.0, 17, 1), -std::sin(angle), std::cos(angle), 151.8, 11)
    };
    const auto rotated_report = rotated_tracker.update(rotated, 0.1);
    if (!rotated_report.updated || rotated_report.ambiguous ||
        rotated_report.pair_subspace_similarity < 0.999) return 7;

    // An exact duplicate candidate creates competing assignments. The tracker
    // must report ambiguity and must not commit the update.
    ModalData ambiguous;
    ambiguous.modes = {mode(4, 101.2, 17, 0), mode(5, 101.2, 17, 0), mode(6, 151.8, 17, 1)};
    const auto uncertain = tracker.update(ambiguous, 0.2);
    if (!uncertain.ambiguous || uncertain.updated) return 8;

    std::cout << "modal_state_tracker_contract: PASS"
              << " score=" << report.best_assignment_score
              << " margin=" << report.assignment_margin
              << " cross_resolution_shape=" << report.tracks[0].shape_similarity
              << " rotated_subspace_similarity=" << rotated_report.pair_subspace_similarity << '\n';
    return 0;
}
