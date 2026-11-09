#include "../../src/fem/cylindrical_bell_model.hpp"
#include "../../src/fem/cylindrical_modal_evaluator.hpp"
#include "../../src/excitation/meridional_influence.hpp"
#include "../../src/excitation/localized_strike.hpp"
#include "round_bell_candidate_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace BellOptimizer;

int main(int argc, char** argv) {
    const std::string output = argc > 1 ? argv[1] : "round_optimizer_integration_probe.csv";
    try {
        // Representative valid round-bell geometry from the current project data.
        BellGeometry geometry;
        geometry.H = 1000.0;
        geometry.R_mouth = 540.7;
        geometry.t_top = 34.7;
        geometry.t_bottom = 70.8;
        geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
        geometry.ellipticity = 1.0;
        if (!geometry.isValid()) throw std::runtime_error("probe geometry is invalid");

        BellMaterial material;
        BellMeshConfig mesh_config;
        mesh_config.nNodes = 15;
        const auto model = FEM::make_cylindrical_bell_model(geometry, material, mesh_config);

        FEM::CylindricalModalEvaluationConfig cfg;
        cfg.harmonics = {0, 1, 2, 3};
        cfg.angular_samples = 64;
        cfg.requested_modes = 18;
        cfg.max_iterations = 1600;
        cfg.tolerance = 1.0e-6;
        cfg.minimum_frequency_hz = 1.0;

        Excitation::MeridionalInfluence legacy_influence;
        ModalData modes;
        if (!FEM::evaluate_cylindrical_bell_modes(
                model, legacy_influence.make(model.mesh_si), cfg, modes))
            throw std::runtime_error("direct cylindrical modal evaluation failed");

        double z_min = model.mesh_si.outer_profile.front().z;
        double z_max = z_min;
        for (const auto& p : model.mesh_si.outer_profile) {
            z_min = std::min(z_min, p.z);
            z_max = std::max(z_max, p.z);
        }
        Excitation::CircularLocalizedStrikeConfig strike;
        strike.center_z = z_min + 0.70 * (z_max - z_min);
        strike.meridional_sigma = 0.012;
        strike.circumferential_sigma = 0.008;
        strike.azimuth_rad = 0.0;
        strike.normal_fraction = 1.0;
        strike.tangent_fraction = 0.0;

        Research::ModalFamilyAdapterConfig family_cfg;
        Research::RoundBellCandidateAnalysis analysis;
        if (!Research::analyzeRoundBellCandidate(
                model.mesh_si, modes, strike, family_cfg, analysis))
            throw std::runtime_error("round-bell candidate analysis rejected FEM/strike data");
        const auto& projections = analysis.signed_strike_projections;
        const auto& families = analysis.modal_families;
        const double max_abs = analysis.maximum_absolute_strike_projection;
        const std::string family_output = output + ".families.csv";
        std::ofstream family_csv(family_output);
        if (!family_csv) throw std::runtime_error("cannot open modal-family CSV");
        family_csv << "family_index,first_mode,mode_count,frequency_hz,dominant_m,dominant_share,combined_strike_projection,harmonic_shares\n";
        family_csv << std::setprecision(12);
        for (std::size_t i = 0; i < families.size(); ++i) {
            const auto& g = families[i];
            family_csv << i << ',' << g.first_mode << ',' << g.mode_count << ','
                       << g.representative_frequency_hz << ',' << g.dominant_azimuthal_order << ','
                       << g.dominant_harmonic_share << ',' << g.combined_strike_projection << ",\"";
            for (std::size_t h = 0; h < g.harmonic_energy_shares.size(); ++h) {
                if (h) family_csv << ';';
                family_csv << g.harmonic_energy_shares[h];
            }
            family_csv << '"' << '\n';
        }
        family_csv.close();

        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("cannot open output CSV");
        csv << "index,frequency_hz,harmonics,projection_signed,projection_abs_normalized,relative_gap_to_next\n";
        csv << std::setprecision(12);
        for (std::size_t i = 0; i < modes.modes.size(); ++i) {
            const auto& m = modes.modes[i];
            double gap = 0.0;
            if (i + 1 < modes.modes.size() && m.frequency_hz > 0.0)
                gap = (modes.modes[i + 1].frequency_hz - m.frequency_hz) / m.frequency_hz;
            csv << i << "," << m.frequency_hz << ",\"";
            for (std::size_t h = 0; h < m.fem_harmonics.size(); ++h) {
                if (h) csv << ';';
                csv << m.fem_harmonics[h];
            }
            csv << "\"," << projections[i] << ',' << std::abs(projections[i]) / max_abs << ',' << gap << '\n';
        }
        csv.close();

        std::cout << "round_optimizer_integration_probe: PASS\n"
                  << "modes=" << modes.modes.size() << "\n"
                  << "strike_center_z_m=" << strike.center_z << "\n"
                  << "meridional_sigma_m=" << strike.meridional_sigma << "\n"
                  << "circumferential_sigma_m=" << strike.circumferential_sigma << "\n"
                  << "max_abs_projection=" << max_abs << "\n"
                  << "modal_family_groups=" << families.size() << "\n"
                  << "family_csv=" << family_output << "\n"
                  << "csv=" << output << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "round_optimizer_integration_probe: FAIL: " << e.what() << '\n';
        return 1;
    }
}
