#pragma once

#include "modal_types.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace BellOptimizer {

class CCXParser {
public:
    static ModalData parse_file(const std::string& filename, bool quiet = true) {
        ModalData data;
        data.source_file = filename;

        std::ifstream file(filename);
        if (!file.is_open()) {
            if (!quiet) std::cerr << "Error: cannot open " << filename << '\n';
            data.source_file = "ERROR: file not found";
            return data;
        }

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(file, line)) lines.push_back(line);

        auto find_section = [&](const char* title) -> size_t {
            for (size_t i = 0; i < lines.size(); ++i)
                if (lines[i].find(title) != std::string::npos) return i;
            return lines.size();
        };

        const size_t eigen_section = find_section("E I G E N V A L U E   O U T P U T");
        const size_t participation_section = find_section("P A R T I C I P A T I O N   F A C T O R S");
        const size_t mass_section = find_section("E F F E C T I V E   M O D A L   M A S S");

        std::vector<ModalMode> candidate_modes;
        double max_abs_eigenvalue = 0.0;
        if (eigen_section != lines.size()) {
            if (!quiet) std::cout << "  Found eigenvalue section\n";
            for (size_t i = eigen_section + 1; i < lines.size(); ++i) {
                const auto& text = lines[i];
                if (text.find("P A R T I C I P A T I O N") != std::string::npos ||
                    text.find("E F F E C T I V E") != std::string::npos) break;
                std::istringstream iss(text);
                int mode_num;
                double eigenvalue, real_part, freq, imag_part;
                if (!(iss >> mode_num >> eigenvalue >> real_part >> freq >> imag_part)) continue;
                if (mode_num < 1) continue;
                ModalMode mode;
                mode.mode_number = mode_num;
                mode.eigenvalue = eigenvalue;
                if (freq > 0.0) mode.frequency_hz = freq;
                else if (real_part > 0.0) mode.frequency_hz = real_part / (2.0 * 3.14159265358979323846);
                if (!(mode.frequency_hz >= 0.0) || !std::isfinite(mode.frequency_hz) ||
                    !std::isfinite(mode.eigenvalue)) continue;
                max_abs_eigenvalue = std::max(max_abs_eigenvalue, std::abs(eigenvalue));
                candidate_modes.push_back(mode);
            }
            const double zero_tolerance = std::max(
                max_abs_eigenvalue * 1e-12, std::numeric_limits<double>::min());
            for (const auto& mode : candidate_modes) {
                // Reject numerical zero/rigid modes relative to the eigenvalue
                // spectrum. Do not discard a physical mode solely because f < 0.5 Hz.
                if (mode.eigenvalue > zero_tolerance && mode.frequency_hz > 0.0)
                    data.modes.push_back(mode);
            }
        }

        auto find_mode = [&](int mode_num) -> ModalMode* {
            for (auto& mode : data.modes)
                if (mode.mode_number == mode_num) return &mode;
            return nullptr;
        };

        if (participation_section != lines.size()) {
            if (!quiet) std::cout << "  Found participation section\n";
            for (size_t i = participation_section + 1; i < lines.size(); ++i) {
                if (lines[i].find("E F F E C T I V E") != std::string::npos) break;
                std::istringstream iss(lines[i]);
                int mode_num;
                double px, py, pz, prx, pry, prz;
                if (!(iss >> mode_num >> px >> py >> pz >> prx >> pry >> prz)) continue;
                if (auto* mode = find_mode(mode_num)) {
                    mode->part_x = px; mode->part_y = py; mode->part_z = pz;
                    mode->part_rx = prx; mode->part_ry = pry; mode->part_rz = prz;
                }
            }
        }

        if (mass_section != lines.size()) {
            if (!quiet) std::cout << "  Found effective modal mass section\n";
            for (size_t i = mass_section + 1; i < lines.size(); ++i) {
                if (lines[i].find("T O T A L") != std::string::npos) break;
                std::istringstream iss(lines[i]);
                int mode_num;
                double mx, my, mz, mrx, mry, mrz;
                if (!(iss >> mode_num >> mx >> my >> mz >> mrx >> mry >> mrz)) continue;
                if (auto* mode = find_mode(mode_num)) {
                    mode->mass_x = mx; mode->mass_y = my; mode->mass_z = mz;
                    mode->mass_rx = mrx; mode->mass_ry = mry; mode->mass_rz = mrz;
                }
            }
        }

        if (eigen_section == lines.size() && participation_section == lines.size() && mass_section == lines.size()) {
            if (!quiet) std::cerr << "Warning: no modal sections found\n";
            data.source_file = "WARNING: no modal data found";
            return data;
        }

        data.sort_by_frequency();
        double max_effective_mass = 0.0;
        for (auto& mode : data.modes) {
            double amp = std::abs(mode.mass_rz);
            if (amp < 1e-10) {
                amp = std::sqrt(mode.mass_x * mode.mass_x + mode.mass_y * mode.mass_y +
                                mode.mass_z * mode.mass_z + mode.mass_rx * mode.mass_rx +
                                mode.mass_ry * mode.mass_ry + mode.mass_rz * mode.mass_rz);
            }
            mode.effective_mass_norm = amp;
            max_effective_mass = std::max(max_effective_mass, amp);
        }
        for (auto& mode : data.modes) {
            const double excitation = std::sqrt(std::max(0.0, mode.effective_mass_norm));
            mode.modal_excitation = excitation;
            mode.acoustic_amplitude = 0.0;
            mode.acoustic_pressure_coefficient = 0.0;
            mode.amplitude = (max_effective_mass > 1e-10)
                ? excitation / std::sqrt(max_effective_mass)
                : 1.0;
            // Legacy CalculiX .dat files do not contain mode shapes, so the
            // surface acoustic coefficient is intentionally unavailable. The
            // mechanical participation remains an explicit compatibility
            // fallback; it must not be interpreted as acoustic pressure.
            mode.acoustic_volume_displacement = 0.0;
        }

        if (!data.modes.empty()) data.fundamental_frequency = data.modes.front().frequency_hz;
        if (!quiet) std::cout << "  Modes: " << data.modes.size() << '\n';
        return data;
    }
};

} // namespace BellOptimizer
