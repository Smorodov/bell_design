#include "excitation/clapper_dynamics.hpp"
#include "core/bell_material.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;
static void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
    try {
        ClapperScalingParams scaling;
        const BellMaterial material = default_material_preset().material;
        const double referenceFrequencies[] = {180.0, 270.0, 410.0};
        const double scales[] = {0.2, 0.5, 1.0, 2.0};
        std::vector<double> fractions, normalizedInitialEnergy;
        for (double scale : scales) {
            ClapperDynamicsConfig config;
            require(makeScaledClapperDynamicsConfig(scale, 0.0, scaling, material, config),
                    "failed to create geometrically scaled clapper");
            const std::vector<Point2D> profile{{0.30 * scale, 0.0}, {0.30 * scale, scale}};
            std::vector<ModalMode> modes(3);
            for (std::size_t i = 0; i < modes.size(); ++i) {
                modes[i].frequency_hz = referenceFrequencies[i] / scale;
                modes[i].modal_mass = 1.0;
                modes[i].inner_surface_harmonics = {0};
                const double radial = 1.0e-3 / std::pow(scale, 1.5);
                modes[i].inner_surface_fourier_radial = {radial, radial};
                modes[i].inner_surface_fourier_axial = {0.0, 0.0};
            }
            ClapperDynamicsResult result;
            require(simulateInternalClapperContact(config, profile, modes, result),
                    "scaled clapper did not complete contact and release");
            require(result.bell_modal_energy_joules > 0.0 && result.bell_energy_transfer_fraction > 0.0,
                    "no positive energy was transferred to bell modes");
            require(result.maximum_relative_energy_error < 2.0e-4,
                    "energy conservation drift exceeded audit tolerance");
            require(std::abs(result.final_energy_joules - result.clapper_energy_at_separation_joules -
                             result.bell_modal_energy_joules) / result.initial_energy_joules < 2.0e-4,
                    "separation energy does not close");
            fractions.push_back(result.bell_energy_transfer_fraction);
            normalizedInitialEnergy.push_back(result.initial_energy_joules / std::pow(scale, 4.0));
            std::cout << "scale=" << scale << " initial_J=" << result.initial_energy_joules
                      << " bell_J=" << result.bell_modal_energy_joules
                      << " transfer=" << result.bell_energy_transfer_fraction
                      << " peak_N=" << result.peak_contact_force_n
                      << " energy_drift=" << result.maximum_relative_energy_error << '\n';
        }
        const auto [minFraction, maxFraction] = std::minmax_element(fractions.begin(), fractions.end());
        require(*minFraction > 0.90 * *maxFraction,
                "energy transfer fraction is not approximately scale invariant");
        const auto [minEnergy, maxEnergy] = std::minmax_element(
            normalizedInitialEnergy.begin(), normalizedInitialEnergy.end());
        require((*maxEnergy - *minEnergy) / *maxEnergy < 1.0e-10,
                "Froude-scaled initial clapper energy does not follow s^4");
        std::cout << "clapper_scale_energy_audit: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "clapper_scale_energy_audit: FAIL: " << e.what() << '\n';
        return 1;
    }
}
