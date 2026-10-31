#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "ccx/ccx_input_writer.hpp"
#include "ccx/ccx_parser.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <regex>

namespace fs = std::filesystem;
namespace {
constexpr int kMeshNodes = 15;
// CalculiX reports a near-zero axial rigid-body mode before the elastic modes.
// Request one extra eigenpair so ten physical modes remain after filtering.
constexpr int kRequestedModes = 11;

bool makeModel(BellMaterial& material, BellGenerator& generator) {
    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};
    generator.setGeometry(geometry);
    const auto validation = generator.validateGeometry();
    if (!validation.passed) {
        std::cerr << "Geometry validation failed: " << validation.reason << '\n';
        return false;
    }
    if (!generator.buildMesh()) {
        std::cerr << "Mesh generation failed.\n";
        return false;
    }
    return true;
}

bool writeInternal(const fs::path& work, const BellMeshData& mesh, const BellMaterial& material) {
    BellOptimizer::InternalModalSolver solver(false, kRequestedModes, true);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{mesh, material, work.string(), "internal_reference", &modal_influence_1};
    if (!solver.solve(problem, data) || data.modes.size() < 5) {
        std::cerr << "Internal FEM failed or returned fewer than 5 physical modes.\n";
        return false;
    }
    std::ofstream csv(work / "internal_modes.csv");
    if (!csv) return false;
    csv << "mode,frequency_hz,mass_rz,amplitude,modal_excitation,acoustic_volume_displacement\n";
    csv << std::setprecision(12);
    std::ofstream shapes(work / "internal_mode_shapes.csv");
    if (!shapes) return false;
    shapes << "mode,node_id,r,z,ur,uz\n" << std::setprecision(12);
    for (size_t i = 0; i < data.modes.size(); ++i) {
        const auto& m = data.modes[i];
        csv << m.mode_number << ',' << m.frequency_hz << ',' << m.mass_rz << ',' << m.amplitude
            << ',' << m.modal_excitation << ',' << m.acoustic_volume_displacement << '\n';
        if (m.nodal_radial.size() != mesh.nodes.size() || m.nodal_axial.size() != mesh.nodes.size()) return false;
        for (size_t j = 0; j < mesh.nodes.size(); ++j) {
            const auto& node = mesh.nodes[j];
            shapes << m.mode_number << ',' << node.id << ',' << node.r << ',' << node.z << ','
                   << m.nodal_radial[j] << ',' << m.nodal_axial[j] << '\n';
        }
    }
    return static_cast<bool>(csv) && static_cast<bool>(shapes);
}

int prepare(const fs::path& work) {
    std::error_code ec;
    fs::create_directories(work, ec);
    if (ec) {
        std::cerr << "Cannot create work directory: " << work << " (" << ec.message() << ")\n";
        return 2;
    }
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = kMeshNodes;
    BellGenerator generator(material, config);
    if (!makeModel(material, generator)) return 3;

    const auto& mesh = generator.getMeshData();
    const fs::path input = work / "bell_current_geometry.inp";
    if (!BellOptimizer::CCXInputWriter::write(mesh, material, input.string(), kRequestedModes)) {
        std::cerr << "Failed to write CalculiX input: " << input << '\n';
        return 4;
    }
    if (!writeInternal(work, mesh, material)) {
        std::cerr << "Failed to write internal FEM reference.\n";
        return 5;
    }
    std::ofstream meta(work / "mesh_info.txt");
    meta << "geometry=corrected_inner_normal_offset\n"
         << "profile_parameter_nodes=" << kMeshNodes << '\n'
         << "mesh_nodes=" << mesh.nodes.size() << '\n'
         << "mesh_elements=" << mesh.elements.size() << '\n'
         << "requested_eigenpairs_including_rigid_modes=" << kRequestedModes << '\n'
         << "material_E=" << std::setprecision(12) << material.E << '\n'
         << "material_nu=" << material.nu << '\n'
         << "material_rho=" << material.rho << '\n';
    if (!meta) return 6;

    std::cout << "Prepared current-geometry CalculiX job.\n"
              << "Mesh nodes: " << mesh.nodes.size() << ", elements: " << mesh.elements.size() << '\n'
              << "Input: " << input.string() << '\n'
              << "Internal reference: " << (work / "internal_modes.csv").string() << '\n';
    return 0;
}

int compare(const fs::path& work) {
    const fs::path dat = work / "bell_current_geometry.dat";
    const auto ccx = BellOptimizer::CCXParser::parse_file(dat.string(), false);
    if (ccx.modes.size() < 10) {
        std::cerr << "CalculiX .dat has fewer than 10 parsed physical modes: " << dat << '\n'
                  << "Check the .dat/.sta/.msg files for solver errors.\n";
        return 7;
    }
    std::ifstream in(work / "internal_modes.csv");
    if (!in) {
        std::cerr << "Missing internal_modes.csv; run prepare stage first.\n";
        return 8;
    }
    std::string line;
    std::getline(in, line); // header
    struct Row { double f = 0.0, mass = 0.0, amp = 0.0; };
    std::vector<Row> internal;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string field;
        Row r;
        if (!std::getline(iss, field, ',')) continue;
        if (!std::getline(iss, field, ',')) continue;
        try { r.f = std::stod(field); } catch (...) { continue; }
        if (!std::getline(iss, field, ',')) continue;
        try { r.mass = std::stod(field); } catch (...) { r.mass = 0.0; }
        if (!std::getline(iss, field, ',')) continue;
        try { r.amp = std::stod(field); } catch (...) { r.amp = 0.0; }
        internal.push_back(r);
    }
    const size_t n = std::min(internal.size(), ccx.modes.size());
    if (n < 10) {
        std::cerr << "Need at least 10 matched physical modes; found " << n << '\n';
        return 9;
    }

    std::ofstream csv(work / "comparison.csv");
    if (!csv) return 10;
    constexpr double kFullRevolutionSectorFactor = 180.0; // CCX axisymmetric output uses a 2-degree sector.
    csv << "mode,internal_hz,ccx_hz,relative_frequency_error,internal_mass_rz,ccx_mass_rz_sector,ccx_mass_rz_full_equivalent,relative_mass_rz_error,internal_amplitude,ccx_amplitude\n";
    csv << std::setprecision(12);
    double max_rel = 0.0;
    double max_mass_rel = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double fi = internal[i].f;
        const double fc = ccx.modes[i].frequency_hz;
        if (!(fi > 0.0) || !(fc > 0.0) || !std::isfinite(fi) || !std::isfinite(fc)) return 11;
        const double rel = std::abs(fi - fc) / fc;
        max_rel = std::max(max_rel, rel);
        const double ccx_mass_full = ccx.modes[i].mass_rz * kFullRevolutionSectorFactor;
        const double mass_rel = (std::abs(internal[i].mass) > 1e-20)
            ? std::abs(internal[i].mass - ccx_mass_full) / std::abs(internal[i].mass) : 0.0;
        max_mass_rel = std::max(max_mass_rel, mass_rel);
        csv << i + 1 << ',' << fi << ',' << fc << ',' << rel << ',' << internal[i].mass << ','
            << ccx.modes[i].mass_rz << ',' << ccx_mass_full << ',' << mass_rel << ','
            << internal[i].amp << ',' << ccx.modes[i].amplitude << '\n';
        std::cout << "mode " << (i + 1) << ": internal=" << fi << " Hz, CCX=" << fc
                  << " Hz, rel.err=" << rel * 100.0 << "%\n";
    }
    // Parse CalculiX nodal displacement eigenvectors from the FRD result. Mode numbers
    // are read from 1PMODE records; node IDs are explicit and need not be assumed ordered.
    struct Shape { double ur = 0.0, uz = 0.0; };
    std::map<int, std::map<int, Shape>> ccx_shapes;
    std::ifstream frd(work / "bell_current_geometry.frd");
    if (!frd) {
        std::cerr << "Missing bell_current_geometry.frd; cannot compare mode shapes.\n";
        return 12;
    }
    std::vector<std::string> frd_lines;
    while (std::getline(frd, line)) frd_lines.push_back(line);
    for (size_t li = 0; li < frd_lines.size(); ++li) {
        if (frd_lines[li].find("1PMODE") == std::string::npos) continue;
        std::istringstream mode_stream(frd_lines[li]);
        std::string token;
        int mode_num = 0;
        while (mode_stream >> token) {
            try { mode_num = std::stoi(token); } catch (...) {}
        }
        if (mode_num <= 0) continue;
        for (size_t j = li + 1; j < frd_lines.size(); ++j) {
            if (frd_lines[j].find("1PMODE") != std::string::npos) break;
            if (frd_lines[j].rfind(" -1", 0) != 0) continue;
            static const std::regex number_re(R"([-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[Ee][-+]?[0-9]+)?)");
            std::vector<double> values;
            for (std::sregex_iterator it(frd_lines[j].begin(), frd_lines[j].end(), number_re), end; it != end; ++it) {
                try { values.push_back(std::stod(it->str())); } catch (...) {}
            }
            if (values.size() < 5) continue;
            const int node_id = static_cast<int>(values[1]);
            ccx_shapes[mode_num][node_id] = {values[2], values[3]};
        }
    }

    std::ifstream internal_shape_file(work / "internal_mode_shapes.csv");
    if (!internal_shape_file) {
        std::cerr << "Missing internal_mode_shapes.csv; rerun prepare stage.\n";
        return 13;
    }
    std::map<int, std::map<int, Shape>> internal_shapes;
    std::getline(internal_shape_file, line);
    while (std::getline(internal_shape_file, line)) {
        std::istringstream row(line);
        std::string field;
        int mode_num = 0, node_id = 0;
        double r = 0.0, z = 0.0, ur = 0.0, uz = 0.0;
        if (!std::getline(row, field, ',')) continue; try { mode_num = std::stoi(field); } catch (...) { continue; }
        if (!std::getline(row, field, ',')) continue; try { node_id = std::stoi(field); } catch (...) { continue; }
        if (!std::getline(row, field, ',')) continue; try { r = std::stod(field); } catch (...) { continue; }
        if (!std::getline(row, field, ',')) continue; try { z = std::stod(field); } catch (...) { continue; }
        if (!std::getline(row, field, ',')) continue; try { ur = std::stod(field); } catch (...) { continue; }
        if (!std::getline(row, field, ',')) continue; try { uz = std::stod(field); } catch (...) { continue; }
        (void)r; (void)z;
        internal_shapes[mode_num][node_id] = {ur, uz};
    }
    std::ofstream shape_csv(work / "mode_shape_comparison.csv");
    if (!shape_csv) return 14;
    shape_csv << "mode,nodes_compared,absolute_correlation,normalized_l2_error,best_fit_scale\n";
    shape_csv << std::setprecision(12);
    double min_corr = 1.0, max_shape_error = 0.0;
    size_t compared_shapes = 0;
    for (const auto& [mode_num, in_nodes] : internal_shapes) {
        const auto ccx_it = ccx_shapes.find(mode_num);
        if (ccx_it == ccx_shapes.end()) continue;
        double dot = 0.0, norm_i = 0.0, norm_c = 0.0;
        size_t count = 0;
        for (const auto& [node_id, a] : in_nodes) {
            const auto b_it = ccx_it->second.find(node_id);
            if (b_it == ccx_it->second.end()) continue;
            const auto& b = b_it->second;
            dot += a.ur * b.ur + a.uz * b.uz;
            norm_i += a.ur * a.ur + a.uz * a.uz;
            norm_c += b.ur * b.ur + b.uz * b.uz;
            ++count;
        }
        if (count == 0 || norm_i <= 0.0 || norm_c <= 0.0) continue;
        const double corr = std::abs(dot) / std::sqrt(norm_i * norm_c);
        const double scale = dot / norm_c; // handles arbitrary sign and magnitude of eigenvectors
        double residual = 0.0;
        for (const auto& [node_id, a] : in_nodes) {
            const auto b_it = ccx_it->second.find(node_id);
            if (b_it == ccx_it->second.end()) continue;
            const auto& b = b_it->second;
            const double dr = a.ur - scale * b.ur;
            const double dz = a.uz - scale * b.uz;
            residual += dr * dr + dz * dz;
        }
        const double l2 = std::sqrt(residual / norm_i);
        shape_csv << mode_num << ',' << count << ',' << corr << ',' << l2 << ',' << scale << '\n';
        min_corr = std::min(min_corr, corr);
        max_shape_error = std::max(max_shape_error, l2);
        ++compared_shapes;
    }
    if (compared_shapes < 10) {
        std::cerr << "Only " << compared_shapes << " mode shapes could be compared; expected 10.\n";
        return 15;
    }

    std::ofstream summary(work / "comparison_summary.txt");
    summary << std::setprecision(12)
            << "compared_modes=" << n << '\n'
            << "max_relative_frequency_error=" << max_rel << '\n'
            << "max_relative_frequency_error_percent=" << max_rel * 100.0 << '\n'
            << "sector_angle_degrees=2\n"
            << "full_revolution_sector_factor=" << kFullRevolutionSectorFactor << '\n'
            << "max_relative_mass_rz_error_after_sector_scaling=" << max_mass_rel << '\n'
            << "compared_mode_shapes=" << compared_shapes << '\n'
            << "minimum_absolute_mode_shape_correlation=" << min_corr << '\n'
            << "maximum_normalized_mode_shape_l2_error=" << max_shape_error << '\n';
    std::cout << "Compared physical modes: " << n << '\n'
              << "Maximum relative frequency error: " << max_rel * 100.0 << "%\n"
              << "Maximum relative mass_rz error after sector scaling: " << max_mass_rel * 100.0 << "%\n"
              << "Compared mode shapes: " << compared_shapes << '\n'
              << "Minimum absolute shape correlation: " << min_corr << '\n'
              << "Maximum normalized shape L2 error: " << max_shape_error << '\n'
              << "CSV: " << (work / "comparison.csv") << '\n'
              << "Shapes: " << (work / "mode_shape_comparison.csv") << '\n'
              << "Summary: " << (work / "comparison_summary.txt") << '\n';
    return 0;
}
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: ccx_current_geometry prepare|compare <work-directory>\n";
        return 1;
    }
    const std::string command = argv[1];
    const fs::path work = fs::absolute(argv[2]);
    if (command == "prepare") return prepare(work);
    if (command == "compare") return compare(work);
    std::cerr << "Unknown command: " << command << '\n';
    return 1;
}
