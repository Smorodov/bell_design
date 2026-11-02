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


static ModalMode full_fourier_mode(int id, double frequency,
                                  const std::vector<int>& harmonics,
                                  int family, bool added_harmonics = false) {
    ModalMode m;
    m.mode_number = id;
    m.frequency_hz = frequency;
    const std::size_t nodes = 4;
    const std::size_t basis_count = harmonics.size();
    m.fem_harmonics = harmonics;
    m.outer_surface_harmonics = harmonics;
    m.fem_node_r = {0.10, 0.20, 0.30, 0.40};
    m.fem_node_z = {0.00, 0.10, 0.20, 0.30};
    m.fem_displacement.assign(nodes * 3u * basis_count, 0.0);
    m.fem_mass_lumped_weights.assign(nodes * 3u * basis_count, 1.0);
    m.outer_surface_fourier_radial.assign(nodes * basis_count, 0.0);
    m.outer_surface_fourier_azimuthal.assign(nodes * basis_count, 0.0);
    m.outer_surface_fourier_axial.assign(nodes * basis_count, 0.0);
    m.outer_surface_radial.assign(nodes, 0.0);
    m.outer_surface_axial.assign(nodes, 0.0);
    for (std::size_t node = 0; node < nodes; ++node) {
        for (std::size_t basis = 0; basis < basis_count; ++basis) {
            double values[3]{};
            if (!(added_harmonics && basis >= 3u)) {
                for (std::size_t component = 0; component < 3u; ++component) {
                    const double phase = 0.31 * static_cast<double>((node + 1u) * (family + 1)) +
                                         0.17 * static_cast<double>((basis + 1u) * (component + 1u));
                    values[component] = std::sin(phase) +
                        0.07 * static_cast<double>(family + 1) * (component + 1u);
                    m.fem_displacement[(node * 3u + component) * basis_count + basis] =
                        values[component];
                }
            }
            const std::size_t index = node * basis_count + basis;
            m.outer_surface_fourier_radial[index] = values[0];
            m.outer_surface_fourier_azimuthal[index] = values[1];
            m.outer_surface_fourier_axial[index] = values[2];
        }
        m.outer_surface_radial[node] = m.outer_surface_fourier_radial[node * basis_count];
        m.outer_surface_axial[node] = m.outer_surface_fourier_axial[node * basis_count];
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

    // Full-volume direct FEM fields use [node][component][basis]. Verify
    // that continuation compares common Fourier components correctly when
    // the candidate adds two new harmonics to the basis.
    const std::vector<int> base_harmonics{0, 1, -1};
    const std::vector<int> enriched_harmonics{0, 1, -1, 2, -2};
    ModalData full_reference;
    full_reference.modes = {
        full_fourier_mode(20, 100.0, base_harmonics, 0),
        full_fourier_mode(21, 150.0, base_harmonics, 1)
    };
    ModalStateTracker full_tracker;
    if (!full_tracker.initialize(full_reference, 0, 1, 0.0)) return 9;
    ModalData full_candidates;
    full_candidates.modes = {
        full_fourier_mode(31, 151.0, enriched_harmonics, 1, true),
        full_fourier_mode(30, 101.0, enriched_harmonics, 0, true)
    };
    const auto full_report = full_tracker.update(full_candidates, 0.1);
    if (!full_report.updated || full_report.ambiguous ||
        full_report.tracks[0].candidate_index != 1 ||
        full_report.tracks[1].candidate_index != 0 ||
        full_report.tracks[0].shape_similarity < 0.999 ||
        full_report.tracks[1].shape_similarity < 0.999) return 10;

    std::cout << "modal_state_tracker_contract: PASS"
              << " score=" << report.best_assignment_score
              << " margin=" << report.assignment_margin
              << " cross_resolution_shape=" << report.tracks[0].shape_similarity
              << " rotated_subspace_similarity=" << rotated_report.pair_subspace_similarity
              << " full_fourier_similarity=" << full_report.tracks[0].shape_similarity << '\n';
    return 0;
}
