#include "../research/round_bell_audit/modal_family_adapter.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::Research;

static ModalMode mode(double f, const std::vector<double>& u) {
    ModalMode m;
    m.frequency_hz = f;
    m.fem_node_r = {0.2};
    m.fem_node_z = {0.4};
    m.fem_harmonics = {0, 1, -1};
    m.fem_displacement = u;
    m.fem_mass_lumped_weights.assign(9, 1.0);
    return m;
}

int main() {
    ModalData data;
    // One cosine and one sine member of the same circular m=1 eigenspace.
    data.modes.push_back(mode(100.0, {0, 1, 0, 0, 0, 0, 0, 0, 0}));
    data.modes.push_back(mode(100.0002, {0, 0, 1, 0, 0, 0, 0, 0, 0}));
    data.modes.push_back(mode(180.0, {1, 0, 0, 0, 0, 0, 0, 0, 0}));
    ModalFamilyAdapterConfig cfg;
    cfg.degeneracy_relative_tolerance = 1.0e-5;
    std::vector<ModalFamilyGroup> groups;
    const std::vector<double> p{3.0, 4.0, 2.0};
    if (!buildModalFamilyGroups(data, p, cfg, groups) || groups.size() != 2) return 1;
    if (groups[0].mode_count != 2 || groups[0].dominant_azimuthal_order != 1 ||
        std::abs(groups[0].combined_strike_projection - 5.0) > 1.0e-12 ||
        std::abs(groups[0].dominant_harmonic_share - 1.0) > 1.0e-12) return 2;

    // Rotate the basis of the degenerate pair by an orthogonal 2x2 matrix.
    const double c = 0.6, s = 0.8;
    ModalData rotated;
    rotated.modes.push_back(mode(100.0, {0, c, s, 0, 0, 0, 0, 0, 0}));
    rotated.modes.push_back(mode(100.0002, {0, -s, c, 0, 0, 0, 0, 0, 0}));
    rotated.modes.push_back(data.modes[2]);
    const std::vector<double> pr{c * 3.0 + s * 4.0, -s * 3.0 + c * 4.0, 2.0};
    std::vector<ModalFamilyGroup> rotated_groups;
    if (!buildModalFamilyGroups(rotated, pr, cfg, rotated_groups) || rotated_groups.size() != 2) return 3;
    if (std::abs(rotated_groups[0].combined_strike_projection - groups[0].combined_strike_projection) > 1.0e-12 ||
        std::abs(rotated_groups[0].dominant_harmonic_share - groups[0].dominant_harmonic_share) > 1.0e-12) return 4;

    // Reject same-sized but geometrically different meshes: modal vectors
    // cannot be compared/grouped safely unless they share the same DOF layout.
    ModalData mismatched = data;
    mismatched.modes[1].fem_node_r[0] += 1.0e-5;
    std::vector<ModalFamilyGroup> rejected;
    if (buildModalFamilyGroups(mismatched, p, cfg, rejected)) return 5;

    // Reject unsorted eigenfrequencies rather than silently grouping bad input.
    ModalData unsorted = data;
    std::swap(unsorted.modes[0], unsorted.modes[2]);
    const std::vector<double> pu{2.0, 4.0, 3.0};
    if (buildModalFamilyGroups(unsorted, pu, cfg, rejected)) return 6;

    std::cout << "modal_family_adapter_contract: PASS\n"
              << "degenerate_pair_projection_norm=" << groups[0].combined_strike_projection << "\n"
              << "orthogonal_rotation_invariance=PASS\n"
              << "separate_family_count=" << groups.size() << "\n"
              << "mesh_identity_guard=PASS\n"
              << "frequency_order_guard=PASS\n";
    return 0;
}
