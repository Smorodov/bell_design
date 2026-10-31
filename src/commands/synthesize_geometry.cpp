#include "synthesize_geometry.hpp"
#include "../core/bell_types.h"
#include "../core/bell_generator.h"
#include "../core/modal_types.hpp"
#include "../core/simple_json.hpp"
#include "../core/math_constants.hpp"
#include "../math/modal_coordinate_morphing.hpp"
#include "../fem/internal_modal_solver.hpp"
#include "../excitation/meridional_influence.hpp"
#include "../synthesis/bell_synthesis.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>

using namespace BellOptimizer;

namespace {
BellGeometry load_geometry(const std::string& filename) {
    const auto data = bell::json::parse_file(filename);
    const auto* g = data.find("geometry");
    if (!g || !g->is_object()) throw std::runtime_error("geometry object is missing");
    BellGeometry out;
    out.H = bell::json::number(*g, "H", 1000.0);
    out.R_mouth = bell::json::number(*g, "R_mouth", 500.0);
    out.t_top = bell::json::number(*g, "t_top", 10.0);
    out.t_bottom = bell::json::number(*g, "t_bottom", 30.0);
    out.ellipticity = bell::json::number(*g, "ellipticity", 1.0);
    const auto* knots = g->find("r_knots");
    if (!knots || !knots->is_array()) throw std::runtime_error("r_knots are missing");
    for (const auto& v : knots->array()) out.r_knots.push_back(v.number(0.0));
    return out;
}

std::vector<Mode> make_modes(const ModalData& data, double& prime) {
    if (data.modes.size() < 5) throw std::runtime_error("at least five modes are required");
    prime = data.modes[1].frequency_hz;
    double max_coeff = 0.0;
    for (const auto& d : data.modes) {
        const double omega = 2.0 * bell::kPi * d.frequency_hz;
        max_coeff = std::max(max_coeff, std::abs(-d.part_rz * d.acoustic_volume_displacement * omega));
    }
    if (max_coeff < 1e-14) max_coeff = 1.0;
    std::vector<Mode> modes;
    for (const auto& d : data.modes) {
        Mode m;
        m.frequency = d.frequency_hz;
        const double omega = 2.0 * bell::kPi * m.frequency;
        const double coeff = -d.part_rz * d.acoustic_volume_displacement * omega;
        m.amplitude = std::abs(coeff) / max_coeff;
        m.phase = coeff < 0.0 ? bell::kPi : 0.0;
        m.decay_time = (m.frequency > 1.0)
            ? std::log(1000.0) / (bell::kPi * m.frequency * 0.0005) : 5.0;
        m.decay_time = std::clamp(m.decay_time, 0.1, 30.0);
        const double ratio = std::max(m.frequency / prime, 1e-9);
        m.attack_time = std::clamp(0.002 + 0.005 / std::sqrt(ratio), 0.001, 0.03);
        modes.push_back(m);
    }
    return modes;
}
}

int command_synthesize_geometry(int argc, char* argv[]) {
    std::string geometry_file = "best_geometry.json";
    std::string output = "bell_geometry.wav";
    std::string workdir;
    double ellipticity = 1.0;
    bool ellipticity_explicit = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            std::cout << "Usage: bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav] [--workdir dir]\n"
                      << "  q=b/a; q=1 is axisymmetric baseline; q<1 compresses the physical y-axis.\n"
                      << "  Frequencies use a reduced-order Rayleigh morph; acoustic modal amplitudes remain baseline estimates.\n";
            return 0;
        } else if (a == "--output" && i + 1 < argc) output = argv[++i];
        else if (a.rfind("--output=", 0) == 0) output = a.substr(9);
        else if (a == "--ellipticity" && i + 1 < argc) { ellipticity = std::stod(argv[++i]); ellipticity_explicit = true; }
        else if (a.rfind("--ellipticity=", 0) == 0) { ellipticity = std::stod(a.substr(14)); ellipticity_explicit = true; }
        else if (a == "--workdir" && i + 1 < argc) workdir = argv[++i];
        else if (a.rfind("--workdir=", 0) == 0) workdir = a.substr(10);
        else if (geometry_file == "best_geometry.json") geometry_file = a;
        else { std::cerr << "Unknown option: " << a << "\n"; return 2; }
    }
    if (!workdir.empty()) {
        std::error_code ec; std::filesystem::current_path(workdir, ec);
        if (ec) { std::cerr << "Cannot change workdir: " << ec.message() << "\n"; return 2; }
    }
    try {
        if (!std::isfinite(ellipticity) || ellipticity <= 0.0 || ellipticity > 1.0)
            throw std::invalid_argument("--ellipticity q must be finite and in (0, 1]");
        const BellGeometry geometry = load_geometry(geometry_file);
        if (!ellipticity_explicit) ellipticity = geometry.ellipticity;
        BellMaterial material;
        BellMeshConfig mesh_config; mesh_config.nNodes = 15;
        BellGenerator generator(material, mesh_config);
        generator.setGeometry(geometry);
        GeometryValidationResult validation;
        if (!generator.buildMesh(&validation)) {
            std::cerr << "Geometry rejected: " << validation.reason << "\n"; return 1;
        }
        InternalModalSolver solver(true, 10, true);
        ModalData data;
        Excitation::MeridionalInfluence influence;
        ModalProblem problem{generator.getMeshData(), material, ".", "geometry_synth", &influence};
        if (!solver.solve(problem, data)) { std::cerr << "Modal solve failed\n"; return 1; }
        const Math::CoordinateScale3D inverse_map{1.0, ellipticity, 1.0};
        if (ellipticity < 1.0) {
            for (auto& mode : data.modes) {
                const auto morph = Math::evaluate_mechanical_morphing(
                    generator.getMeshData(), material, mode, inverse_map);
                mode.frequency_hz = morph.morphed_frequency_hz;
                mode.eigenvalue = std::pow(2.0 * bell::kPi * mode.frequency_hz, 2.0);
            }
            data.sort_by_frequency();
        }
        double prime = 0.0;
        auto modes = make_modes(data, prime);
        SynthesisConfig cfg;
        cfg.duration = 8.0;
        cfg.reverb_mix = 0.0;
        cfg.target_fundamental = prime;
        BellSynthesizer synth(cfg);
        GeometryParams geo{};
        geo.height = geometry.H / 1000.0;
        geo.r_max = geometry.R_mouth / 1000.0;
        geo.thickness = geometry.t_bottom / 1000.0;
        synth.generate_modes(geo, modes);
        synth.synthesize();
        if (!synth.save_wav(output)) { std::cerr << "Failed to write " << output << "\n"; return 1; }
        std::vector<double> frequencies; for (const auto& m : modes) frequencies.push_back(m.frequency);
        const std::filesystem::path p(output);
        const auto stem = p.stem().string();
        synth.save_spectrum_csv(stem + "_spectrum.csv", frequencies);
        synth.save_all_envelopes_csv(stem + "_envelopes.csv");
        std::cout << std::fixed << std::setprecision(2)
                  << "Generated " << output << ", q=b/a=" << ellipticity
                  << ", Prime=" << prime << " Hz, duration="
                  << synth.get_actual_duration() << " s\n";
        std::ofstream metadata(std::filesystem::path(output).replace_extension("_morphing.txt"));
        if (metadata) {
            metadata << "BellOptimizer coordinate-morphing listening sample\n"
                     << "q=b/a: " << ellipticity << "\n"
                     << "Frequency model: single-mode Rayleigh quotient on transformed baseline mode shapes\n"
                     << "Acoustic modal amplitudes and phases: baseline estimates, not recalculated for the non-axisymmetric field\n"
                     << "This is a listening preview, not a validated prediction of a physical elliptical bell.\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "synth-geometry error: " << e.what() << "\n"; return 1;
    }
}
