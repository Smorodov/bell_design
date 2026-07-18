// =============================================================================
// ccx_parser.hpp - ПАРСЕР .dat ФАЙЛОВ CALCULIX
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <cctype>

#include <fmt/format.h>

namespace BellOptimizer {

struct CCXModeData {
    int mode_number = 0;
    double eigenvalue = 0.0;
    double frequency_hz = 0.0;
    double amplitude = 0.0;
    double decay_time = 0.0;

    double part_x = 0.0, part_y = 0.0, part_z = 0.0;
    double part_rx = 0.0, part_ry = 0.0, part_rz = 0.0;

    double mass_x = 0.0, mass_y = 0.0, mass_z = 0.0;
    double mass_rx = 0.0, mass_ry = 0.0, mass_rz = 0.0;
    
    double effective_mass_norm = 0.0;
};

struct CCXData {
    std::vector<CCXModeData> modes;
    double fundamental_frequency = 0.0;
    std::string source_file;
    
    double sample_rate = 48000.0;
    double duration = 3.0;
    double t60_base = 2.0;
    double target_rms = 0.70;
    std::string material;

    void sort_by_frequency() {
        std::sort(modes.begin(), modes.end(),
                  [](const CCXModeData& a, const CCXModeData& b) {
                      return a.frequency_hz < b.frequency_hz;
                  });
    }

    void normalize_amplitudes(double max_amp = 1.0) {
        if (modes.empty()) return;
        double max_amp_found = 0.0;
        for (const auto& mode : modes) {
            if (mode.amplitude > max_amp_found) max_amp_found = mode.amplitude;
        }
        if (max_amp_found > 1e-10) {
            double scale = max_amp / max_amp_found;
            for (auto& mode : modes) mode.amplitude *= scale;
        }
    }

    size_t size() const { return modes.size(); }
    bool empty() const { return modes.empty(); }
};

class CCXParser {
public:
    static CCXData parse_file(const std::string& filename, bool quiet = true) {
        CCXData data;
        data.source_file = filename;

        std::ifstream file(filename);
        if (!file.is_open()) {
            if (!quiet) fmt::print(stderr, "❌ Ошибка: не удалось открыть файл {}\n", filename);
            data.source_file = "ERROR: file not found";
            return data;
        }

        if (!quiet) fmt::print("📂 Загрузка данных из {}...\n", filename);

        std::string line;
        bool found_eigenvalue_section = false;
        bool found_participation_section = false;
        bool found_effective_mass_section = false;
        
        while (std::getline(file, line)) {
            if (line.find("E I G E N V A L U E   O U T P U T") != std::string::npos) {
                if (!quiet) fmt::print("   📊 Найдена секция EIGENVALUE OUTPUT\n");
                found_eigenvalue_section = true;
                for (int i = 0; i < 5 && std::getline(file, line); ++i) {}
                
                while (std::getline(file, line)) {
                    if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) break;
                    
                    std::istringstream iss(line);
                    int mode_num;
                    double eigenvalue, real_part, freq, imag_part;
                    
                    if (iss >> mode_num >> eigenvalue >> real_part >> freq >> imag_part) {
                        CCXModeData mode;
                        mode.mode_number = mode_num;
                        mode.eigenvalue = eigenvalue;
                        
                        if (freq > 0) {
                            mode.frequency_hz = freq;
                        } else if (real_part > 0) {
                            mode.frequency_hz = real_part / (2.0 * M_PI);
                        }
                        
                        if (mode.frequency_hz > 0.5) {
                            data.modes.push_back(mode);
                        }
                    }
                }
            }
            
            if (line.find("P A R T I C I P A T I O N   F A C T O R S") != std::string::npos) {
                if (!quiet) fmt::print("   📊 Найдена секция PARTICIPATION FACTORS\n");
                found_participation_section = true;
                for (int i = 0; i < 3 && std::getline(file, line); ++i) {}
                
                while (std::getline(file, line)) {
                    if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) break;
                    
                    std::istringstream iss(line);
                    int mode_num;
                    double px, py, pz, prx, pry, prz;
                    
                    if (iss >> mode_num >> px >> py >> pz >> prx >> pry >> prz) {
                        for (auto& mode : data.modes) {
                            if (mode.mode_number == mode_num) {
                                mode.part_x = px;
                                mode.part_y = py;
                                mode.part_z = pz;
                                mode.part_rx = prx;
                                mode.part_ry = pry;
                                mode.part_rz = prz;
                                break;
                            }
                        }
                    }
                }
            }
            
            if (line.find("E F F E C T I V E   M O D A L   M A S S") != std::string::npos) {
                if (!quiet) fmt::print("   📊 Найдена секция EFFECTIVE MODAL MASS\n");
                found_effective_mass_section = true;
                for (int i = 0; i < 3 && std::getline(file, line); ++i) {}
                
                while (std::getline(file, line)) {
                    if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) break;
                    
                    std::istringstream iss(line);
                    int mode_num;
                    double mx, my, mz, mrx, mry, mrz;
                    
                    if (iss >> mode_num >> mx >> my >> mz >> mrx >> mry >> mrz) {
                        for (auto& mode : data.modes) {
                            if (mode.mode_number == mode_num) {
                                mode.mass_x = mx;
                                mode.mass_y = my;
                                mode.mass_z = mz;
                                mode.mass_rx = mrx;
                                mode.mass_ry = mry;
                                mode.mass_rz = mrz;
                                break;
                            }
                        }
                    }
                }
            }
        }
        
        file.close();

        // Проверяем, что файл содержит данные
        if (!found_eigenvalue_section && !found_participation_section && !found_effective_mass_section) {
            if (!quiet) fmt::print(stderr, "⚠️ Предупреждение: файл не содержит секций с модальными данными\n");
            data.source_file = "WARNING: no modal data found";
            return data;
        }

        data.sort_by_frequency();

        double max_amp = 0.0;
        for (auto& mode : data.modes) {
            double amp = std::abs(mode.mass_rz);
            if (amp < 1e-10) {
                amp = std::sqrt(
                    mode.mass_x * mode.mass_x +
                    mode.mass_y * mode.mass_y +
                    mode.mass_z * mode.mass_z +
                    mode.mass_rx * mode.mass_rx +
                    mode.mass_ry * mode.mass_ry +
                    mode.mass_rz * mode.mass_rz
                );
            }
            mode.effective_mass_norm = amp;
            if (amp > max_amp) max_amp = amp;
        }
        
        if (max_amp > 1e-10) {
            for (auto& mode : data.modes) {
                mode.amplitude = mode.effective_mass_norm / max_amp;
            }
        } else {
            for (auto& mode : data.modes) {
                mode.amplitude = 1.0;
            }
        }

        if (!data.modes.empty()) {
            data.fundamental_frequency = data.modes[0].frequency_hz;
        }

        if (!quiet) fmt::print("   📊 Оставлено мод: {}\n", data.modes.size());
        return data;
    }
};

} // namespace BellOptimizer