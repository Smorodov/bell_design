#include "modal_family_adapter.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace BellOptimizer::Research {
namespace {

bool validMode(const ModalMode& mode, std::size_t& nodes, std::size_t& basis) {
    nodes = mode.fem_node_r.size();
    basis = mode.fem_harmonics.size();
    if (nodes == 0 || basis == 0 || mode.fem_node_z.size() != nodes ||
        mode.fem_displacement.size() != nodes * 3u * basis ||
        mode.fem_mass_lumped_weights.size() != mode.fem_displacement.size() ||
        !std::isfinite(mode.frequency_hz) || mode.frequency_hz <= 0.0) return false;
    for (double v : mode.fem_displacement) if (!std::isfinite(v)) return false;
    for (double v : mode.fem_mass_lumped_weights) if (!std::isfinite(v) || v < 0.0) return false;
    return true;
}

bool harmonicShares(const ModalMode& mode, std::vector<double>& shares) {
    std::size_t nodes = 0, basis = 0;
    if (!validMode(mode, nodes, basis)) return false;
    std::size_t max_order = 0;
    for (int h : mode.fem_harmonics) {
        const auto order = static_cast<std::size_t>(h < 0 ? -static_cast<long long>(h) : h);
        max_order = std::max(max_order, order);
    }
    if (max_order > 4096u) return false; // reject malformed/unreasonably large bases
    std::vector<double> energy(max_order + 1u, 0.0);
    for (std::size_t n = 0; n < nodes; ++n) {
        for (std::size_t c = 0; c < 3; ++c) {
            for (std::size_t b = 0; b < basis; ++b) {
                const std::size_t i = (n * 3u + c) * basis + b;
                const double u = mode.fem_displacement[i];
                const int harmonic = mode.fem_harmonics[b];
                const auto order = static_cast<std::size_t>(harmonic < 0
                    ? -static_cast<long long>(harmonic) : harmonic);
                energy[order] += mode.fem_mass_lumped_weights[i] * u * u;
            }
        }
    }
    double total = 0.0;
    for (double e : energy) total += e;
    if (!(total > std::numeric_limits<double>::min()) || !std::isfinite(total)) return false;
    shares.resize(energy.size());
    for (std::size_t i = 0; i < energy.size(); ++i) shares[i] = energy[i] / total;
    return true;
}

} // namespace

bool buildModalFamilyGroups(const ModalData& data,
                            const std::vector<double>& projections,
                            const ModalFamilyAdapterConfig& cfg,
                            std::vector<ModalFamilyGroup>& groups) {
    groups.clear();
    const std::size_t count = data.modes.size();
    if (count == 0 || projections.size() != count ||
        !std::isfinite(cfg.degeneracy_relative_tolerance) || cfg.degeneracy_relative_tolerance < 0.0 ||
        !std::isfinite(cfg.minimum_harmonic_share) || cfg.minimum_harmonic_share < 0.0 ||
        cfg.minimum_harmonic_share > 1.0) return false;

    std::vector<std::vector<double>> shares(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (!std::isfinite(projections[i]) || !harmonicShares(data.modes[i], shares[i])) return false;
        if (i && data.modes[i].frequency_hz < data.modes[i - 1].frequency_hz) return false;
        if (i) {
            const auto& reference = data.modes[0];
            const auto& current = data.modes[i];
            if (current.fem_harmonics != reference.fem_harmonics ||
                current.fem_node_r.size() != reference.fem_node_r.size() ||
                current.fem_node_z.size() != reference.fem_node_z.size()) return false;
            for (std::size_t n = 0; n < reference.fem_node_r.size(); ++n) {
                if (!std::isfinite(reference.fem_node_r[n]) || !std::isfinite(reference.fem_node_z[n]) ||
                    !std::isfinite(current.fem_node_r[n]) || !std::isfinite(current.fem_node_z[n]) ||
                    std::abs(current.fem_node_r[n] - reference.fem_node_r[n]) > 1.0e-12 ||
                    std::abs(current.fem_node_z[n] - reference.fem_node_z[n]) > 1.0e-12) return false;
            }
        }
    }

    std::size_t first = 0;
    while (first < count) {
        std::size_t last = first + 1;
        while (last < count) {
            const double f0 = data.modes[last - 1].frequency_hz;
            const double f1 = data.modes[last].frequency_hz;
            const double relative_gap = (f1 - f0) / std::max(f0, 1.0);
            if (relative_gap > cfg.degeneracy_relative_tolerance) break;
            ++last;
        }

        ModalFamilyGroup group;
        group.first_mode = first;
        group.mode_count = last - first;
        group.representative_frequency_hz = 0.0;
        group.harmonic_energy_shares.assign(shares[first].size(), 0.0);
        double projection_sq = 0.0;
        for (std::size_t i = first; i < last; ++i) {
            group.representative_frequency_hz += data.modes[i].frequency_hz;
            projection_sq += projections[i] * projections[i];
            for (std::size_t h = 0; h < shares[i].size(); ++h)
                group.harmonic_energy_shares[h] += shares[i][h];
        }
        const double inv_n = 1.0 / static_cast<double>(group.mode_count);
        group.representative_frequency_hz *= inv_n;
        for (double& e : group.harmonic_energy_shares) e *= inv_n;
        group.projection_energy = projection_sq;
        group.combined_strike_projection = std::sqrt(projection_sq);
        const auto dominant = std::max_element(group.harmonic_energy_shares.begin(),
                                               group.harmonic_energy_shares.end());
        group.dominant_azimuthal_order = static_cast<int>(dominant - group.harmonic_energy_shares.begin());
        group.dominant_harmonic_share = *dominant;
        if (group.dominant_harmonic_share < cfg.minimum_harmonic_share)
            group.dominant_azimuthal_order = -1; // mixed/ambiguous family; do not force a label
        groups.push_back(std::move(group));
        first = last;
    }
    return true;
}

} // namespace BellOptimizer::Research
