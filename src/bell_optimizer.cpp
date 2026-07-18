// bell_optimizer.cpp
#include "bell_optimizer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>

#include <fmt/format.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// ============================================================
// СТРУКТУРА ДЛЯ СОХРАНЕНИЯ РЕЗУЛЬТАТА
// ============================================================
struct BestResult {
    double error = 0.0;
    BellGeometry geometry;
    CCXData ccx_data;
    
    bool save(const std::string& filename) const {
        json data;
        
        data["error"] = error;
        
        data["geometry"]["H"] = geometry.H;
        data["geometry"]["R_mouth"] = geometry.R_mouth;
        data["geometry"]["t_top"] = geometry.t_top;
        data["geometry"]["t_bottom"] = geometry.t_bottom;
        data["geometry"]["r_knots"] = geometry.r_knots;
        
        if (ccx_data.modes.size() > 0) {
            data["frequencies"]["hum"] = ccx_data.modes[0].frequency_hz;
        }
        if (ccx_data.modes.size() > 1) {
            data["frequencies"]["prime"] = ccx_data.modes[1].frequency_hz;
        }
        if (ccx_data.modes.size() > 2) {
            data["frequencies"]["tertia"] = ccx_data.modes[2].frequency_hz;
        }
        if (ccx_data.modes.size() > 3) {
            data["frequencies"]["quinta"] = ccx_data.modes[3].frequency_hz;
        }
        if (ccx_data.modes.size() > 4) {
            data["frequencies"]["nominal"] = ccx_data.modes[4].frequency_hz;
        }
        
        std::ofstream file(filename);
        if (!file.is_open()) return false;
        
        file << data.dump(4);
        file.close();
        return true;
    }
    
    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) return false;
        
        try {
            json data = json::parse(file);
            
            error = data.value("error", 0.0);
            
            if (data.contains("geometry")) {
                auto& g = data["geometry"];
                geometry.H = g.value("H", 1000.0);
                geometry.R_mouth = g.value("R_mouth", 400.0);
                geometry.t_top = g.value("t_top", 10.0);
                geometry.t_bottom = g.value("t_bottom", 60.0);
                
                if (g.contains("r_knots") && g["r_knots"].is_array()) {
                    geometry.r_knots.clear();
                    for (auto& item : g["r_knots"]) {
                        geometry.r_knots.push_back(item.get<double>());
                    }
                }
            }
            
            file.close();
            return true;
        } catch (const std::exception& e) {
            fmt::print(stderr, "❌ Ошибка загрузки JSON: {}\n", e.what());
            file.close();
            return false;
        }
    }
};

// ============================================================
// КОНСТРУКТОР
// ============================================================
EvolutionaryOptimizer::EvolutionaryOptimizer(const BellMaterial& material, 
                                             const TargetData& targets,
                                             const EvolutionParams& evo_params)
    : material_(material)
    , targets_(targets)
    , evo_params_(evo_params)
    , rng_(rd_())
{
}

// ============================================================
// СТАНДАРТНАЯ ГЕОМЕТРИЯ (РЕЗЕРВНЫЙ ВАРИАНТ)
// ============================================================
BellGeometry EvolutionaryOptimizer::getDefaultGeometry() {
    BellGeometry geom;
    geom.H = 1000.0;
    geom.R_mouth = 500.0;
    geom.t_top = 10.0;
    geom.t_bottom = 30.0;
    geom.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    return geom;
}

// ============================================================
// РАЗУМНЫЕ НАЧАЛЬНЫЕ ГЕОМЕТРИИ
// ============================================================
std::vector<BellGeometry> EvolutionaryOptimizer::getReasonableGeometries() {
    std::vector<BellGeometry> geoms;
    
    BellGeometry g1;
    g1.H = 1000.0;
    g1.R_mouth = 500.0;
    g1.t_top = 10.0;
    g1.t_bottom = 30.0;
    g1.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    geoms.push_back(g1);
    
    BellGeometry g2;
    g2.H = 1000.0;
    g2.R_mouth = 600.0;
    g2.t_top = 12.0;
    g2.t_bottom = 35.0;
    g2.r_knots = {0.0, 0.12, 0.35, 0.65, 1.0};
    geoms.push_back(g2);
    
    BellGeometry g3;
    g3.H = 1000.0;
    g3.R_mouth = 400.0;
    g3.t_top = 8.0;
    g3.t_bottom = 25.0;
    g3.r_knots = {0.0, 0.18, 0.45, 0.75, 1.0};
    geoms.push_back(g3);
    
    BellGeometry g4;
    g4.H = 1500.0;
    g4.R_mouth = 750.0;
    g4.t_top = 15.0;
    g4.t_bottom = 45.0;
    g4.r_knots = {0.0, 0.14, 0.38, 0.68, 1.0};
    geoms.push_back(g4);
    
    BellGeometry g5;
    g5.H = 600.0;
    g5.R_mouth = 300.0;
    g5.t_top = 6.0;
    g5.t_bottom = 18.0;
    g5.r_knots = {0.0, 0.16, 0.42, 0.72, 1.0};
    geoms.push_back(g5);
    
    return geoms;
}

// ============================================================
// УСТАНОВИТЬ НАЧАЛЬНОЕ ПРИБЛИЖЕНИЕ
// ============================================================
void EvolutionaryOptimizer::setInitialGuess(const BellGeometry& initial) {
    initial_guess_ = initial;
    has_initial_guess_ = true;
}

// ============================================================
// ЗАГРУЗКА ИЗ ФАЙЛА С ПЕРЕСЧЁТОМ ПОД НОВУЮ ЦЕЛЕВУЮ ЧАСТОТУ
// ============================================================
bool EvolutionaryOptimizer::loadBestFromFile(const std::string& filename) {
    BestResult result;
    if (!result.load(filename)) {
        return false;
    }
    
    BellGeometry loaded_geometry = result.geometry;
    
    fmt::print("📂 Loaded geometry from {}\n", filename);
    fmt::print("   Old H = {:.1f}, R_mouth = {:.1f}\n", 
               loaded_geometry.H, loaded_geometry.R_mouth);
    fmt::print("   Old error = {:.6f}\n", result.error);
    fmt::print("   🔄 Rescaling to target frequency: {:.1f} Hz\n", 
               targets_.prime_freq_target);
    
    // ============================================================
    // 1. РАСЧИТЫВАЕМ ЗАГРУЖЕННУЮ ГЕОМЕТРИЮ
    // ============================================================
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material_, config);
    generator.setGeometry(loaded_geometry);
    
    if (!generator.buildMesh()) {
        fmt::print("   ❌ Failed to build mesh\n");
        return false;
    }
    
    if (!generator.saveToFile("bell.inp")) {
        fmt::print("   ❌ Failed to save bell.inp\n");
        return false;
    }
    
#ifdef _WIN32
    int ccx_result = std::system("ccx bell > nul 2>&1");
#else
    int ccx_result = std::system("ccx bell > /dev/null 2>&1");
#endif
    if (ccx_result != 0) {
        fmt::print("   ❌ CalculiX failed\n");
        return false;
    }
    
    CCXData first_data = BellOptimizer::CCXParser::parse_file("bell.dat", true);
    if (first_data.modes.empty()) {
        fmt::print("   ❌ No modes found\n");
        return false;
    }
    
    double prime_freq = 0.0;
    if (first_data.modes.size() > 1) {
        prime_freq = first_data.modes[1].frequency_hz;
    }
    
    if (prime_freq < 1.0) {
        fmt::print("   ❌ Invalid prime frequency: {:.2f}\n", prime_freq);
        return false;
    }
    
    fmt::print("   Current Prime: {:.2f} Hz\n", prime_freq);
    
    // ============================================================
    // 2. МАСШТАБИРУЕМ ПОД НОВУЮ ЦЕЛЕВУЮ ЧАСТОТУ
    // ============================================================
    double scale = prime_freq / targets_.prime_freq_target;
    
    if (scale < 0.05 || scale > 20.0) {
        fmt::print("   ❌ Invalid scale: {:.3f}\n", scale);
        return false;
    }
    
    fmt::print("   Scale factor: {:.3f}\n", scale);
    
    BellGeometry scaled_geometry = loaded_geometry;
    scaled_geometry.scale(scale);
    
    // ============================================================
    // 3. ПРОВЕРКИ ПОСЛЕ МАСШТАБИРОВАНИЯ (ТОЛЬКО ПРОПОРЦИИ)
    // ============================================================
    
    double t_top_ratio = scaled_geometry.t_top / scaled_geometry.H;
    double t_bottom_ratio = scaled_geometry.t_bottom / scaled_geometry.H;
    
    if (t_top_ratio < 0.005 || t_top_ratio > 0.050) {
        fmt::print("   ❌ t_top/H ratio out of range: {:.4f}\n", t_top_ratio);
        return false;
    }
    if (t_bottom_ratio < 0.010 || t_bottom_ratio > 0.100) {
        fmt::print("   ❌ t_bottom/H ratio out of range: {:.4f}\n", t_bottom_ratio);
        return false;
    }
    
    double diameter = 2.0 * scaled_geometry.R_mouth;
    double ratio_diam = diameter / scaled_geometry.H;
    if (ratio_diam < 0.7 || ratio_diam > 1.3) {
        fmt::print("   ❌ Diameter/Height ratio out of range: {:.3f}\n", ratio_diam);
        return false;
    }
    
    fmt::print("   Scaled H = {:.1f}, R_mouth = {:.1f}\n", scaled_geometry.H, scaled_geometry.R_mouth);
    fmt::print("   t_top/H = {:.4f}, t_bottom/H = {:.4f}\n", t_top_ratio, t_bottom_ratio);
    
    // ============================================================
    // 4. ПЕРЕСЧИТЫВАЕМ МАСШТАБИРОВАННУЮ ГЕОМЕТРИЮ
    // ============================================================
    generator.setGeometry(scaled_geometry);
    
    if (!generator.buildMesh()) {
        fmt::print("   ❌ Failed to build mesh for scaled geometry\n");
        return false;
    }
    
    if (!generator.saveToFile("bell.inp")) {
        fmt::print("   ❌ Failed to save bell.inp for scaled geometry\n");
        return false;
    }
    
#ifdef _WIN32
    ccx_result = std::system("ccx bell > nul 2>&1");
#else
    ccx_result = std::system("ccx bell > /dev/null 2>&1");
#endif
    if (ccx_result != 0) {
        fmt::print("   ❌ CalculiX failed for scaled geometry\n");
        return false;
    }
    
    CCXData scaled_data = BellOptimizer::CCXParser::parse_file("bell.dat", true);
    if (scaled_data.modes.empty()) {
        fmt::print("   ❌ No modes found for scaled geometry\n");
        return false;
    }
    
    // ============================================================
    // 5. СОХРАНЯЕМ ПЕРЕСЧИТАННЫЕ ДАННЫЕ
    // ============================================================
    setInitialGuess(scaled_geometry);
    best_error_ = computeError(scaled_data, scaled_geometry);
    best_geometry_ = scaled_geometry;
    best_ccx_data_ = scaled_data;
    
    fmt::print("✅ Rescaled successfully!\n");
    fmt::print("   New error: {:.6f}\n", best_error_);
    fmt::print("   New Prime: {:.2f} Hz\n", best_ccx_data_.modes[1].frequency_hz);
    
    saveBestResult();
    
    return true;
}

// ============================================================
// КОПИРОВАНИЕ ФАЙЛА
// ============================================================
void EvolutionaryOptimizer::copyFile(const std::string& src, const std::string& dst) {
    std::ifstream src_file(src, std::ios::binary);
    if (!src_file.is_open()) return;
    std::ofstream dst_file(dst, std::ios::binary);
    if (!dst_file.is_open()) return;
    dst_file << src_file.rdbuf();
}

// ============================================================
// СОХРАНЕНИЕ ЛУЧШЕГО РЕЗУЛЬТАТА
// ============================================================
void EvolutionaryOptimizer::saveBestResult() {
    fmt::print("   🔍 Saving best result...\n");
    fmt::print("      Error: {:.6f}\n", best_error_);
    fmt::print("      Modes: {}\n", best_ccx_data_.modes.size());
    
    if (best_error_ >= 1e8) {
        fmt::print("      ❌ Error too high (>= 1e8)\n");
        return;
    }
    
    if (best_ccx_data_.modes.size() < 3) {
        fmt::print("      ❌ Not enough modes: {} < 3\n", best_ccx_data_.modes.size());
        return;
    }
    
    bool has_valid_freq = false;
    for (const auto& mode : best_ccx_data_.modes) {
        if (mode.frequency_hz > 20.0) {
            has_valid_freq = true;
            break;
        }
    }
    
    if (!has_valid_freq) {
        fmt::print("      ❌ No valid frequencies (> 20 Hz)\n");
        return;
    }
    
    BestResult result;
    result.error = best_error_;
    result.geometry = best_geometry_;
    result.ccx_data = best_ccx_data_;
    
    fmt::print("      💾 Writing to best_geometry.json...\n");
    if (result.save("best_geometry.json")) {
        fmt::print("      ✅ Saved successfully!\n");
    } else {
        fmt::print("      ❌ SAVE FAILED! Cannot open file for writing.\n");
        fmt::print("      ⚠️  Check file permissions or if file is open in another program.\n");
    }
    
    copyFile("bell.inp", "best_bell.inp");
    copyFile("bell.dat", "best_bell.dat");
    fmt::print("      ✅ Copied bell.inp -> best_bell.inp, bell.dat -> best_bell.dat\n");
}

// ============================================================
// ГЕНЕРАЦИЯ МОНОТОННЫХ РАДИУСОВ
// ============================================================
std::vector<double> EvolutionaryOptimizer::generateMonotonicRadii(double Rmouth, int num_internal) {
    std::vector<double> radii;
    radii.reserve(num_internal + 2);
    
    radii.push_back(0.0);
    
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<double> fractions(num_internal);
    for (int i = 0; i < num_internal; ++i) {
        fractions[i] = dist(rng_);
    }
    std::sort(fractions.begin(), fractions.end());
    
    for (int i = 0; i < num_internal; ++i) {
        double min_fraction = 0.05 + i * 0.25;
        double max_fraction = 0.15 + i * 0.28;
        double fraction = min_fraction + fractions[i] * (max_fraction - min_fraction);
        radii.push_back(fraction * Rmouth);
    }
    
    radii.push_back(Rmouth);
    return radii;
}

// ============================================================
// ГЕНЕРАЦИЯ СЛУЧАЙНОЙ ГЕОМЕТРИИ
// ============================================================
BellGeometry EvolutionaryOptimizer::generateRandomGeometry() {
    std::uniform_real_distribution<double> distH(500.0, 2000.0);
    std::uniform_real_distribution<double> distR(0.0, 1.0);
    
    BellGeometry geom;
    
    geom.H = distH(rng_);
    
    double r_min = 0.4 * geom.H;
    double r_max = 0.8 * geom.H;
    geom.R_mouth = r_min + distR(rng_) * (r_max - r_min);
    
    std::uniform_real_distribution<double> distTtopRatio(0.008, 0.025);
    std::uniform_real_distribution<double> distTbottomRatio(0.03, 0.08);
    
    geom.t_top = distTtopRatio(rng_) * geom.H;
    geom.t_bottom = distTbottomRatio(rng_) * geom.H;
    
    int num_internal = evo_params_.num_knots - 2;
    std::vector<double> r_abs = generateMonotonicRadii(geom.R_mouth, num_internal);
    
    geom.r_knots.clear();
    for (double r : r_abs) {
        geom.r_knots.push_back(r / geom.R_mouth);
    }
    
    return geom;
}

// ============================================================
// СОЗДАНИЕ СЛУЧАЙНОГО КАНДИДАТА
// ============================================================
EvolutionaryOptimizer::Candidate EvolutionaryOptimizer::createRandomCandidate() {
    Candidate cand;
    cand.geometry = generateRandomGeometry();
    cand.error = 1e9;
    cand.valid = false;
    return cand;
}

// ============================================================
// СОЗДАНИЕ КАНДИДАТА ВБЛИЗИ ЗАГРУЖЕННОЙ ГЕОМЕТРИИ (АДАПТИВНЫЙ РАДИУС)
// ============================================================
EvolutionaryOptimizer::Candidate EvolutionaryOptimizer::createNearCandidate(const BellGeometry& base, double error) {
    Candidate cand;
    cand.error = 1e9;
    cand.valid = false;
    
    double max_error = 1.0;
    double scale = std::min(error / max_error, 1.0);
    double radius = 0.03 + 0.17 * scale;
    
    std::uniform_real_distribution<double> dist(-radius, radius);
    
    const int MAX_ATTEMPTS = 100000;
    
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        BellGeometry geom = base;
        
        double delta_H = dist(rng_) * base.H;
        geom.H = std::max(500.0, std::min(2000.0, base.H + delta_H));
        
        double r_min = 0.4 * geom.H;
        double r_max = 0.8 * geom.H;
        double delta_R = dist(rng_) * base.R_mouth;
        geom.R_mouth = std::max(r_min, std::min(r_max, base.R_mouth + delta_R));
        
        double delta_t_top = dist(rng_) * base.t_top;
        geom.t_top = std::max(2.0, std::min(60.0, base.t_top + delta_t_top));
        
        double delta_t_bottom = dist(rng_) * base.t_bottom;
        geom.t_bottom = std::max(10.0, std::min(160.0, base.t_bottom + delta_t_bottom));
        
        for (size_t i = 1; i < geom.r_knots.size() - 1; ++i) {
            double delta_r = dist(rng_) * 0.15;
            geom.r_knots[i] = std::max(0.0, std::min(1.0, base.r_knots[i] + delta_r));
        }
        
        for (size_t i = 1; i < geom.r_knots.size(); ++i) {
            if (geom.r_knots[i] <= geom.r_knots[i-1]) {
                geom.r_knots[i] = geom.r_knots[i-1] + 0.01;
            }
            if (geom.r_knots[i] > 1.0) geom.r_knots[i] = 1.0;
        }
        if (geom.r_knots.back() != 1.0) geom.r_knots.back() = 1.0;
        
        if (!geom.isValid()) continue;
        
        BellMeshConfig config;
        config.nNodes = 15;
        BellGenerator generator(material_, config);
        generator.setGeometry(geom);
        
        auto validation = generator.validateGeometry();
        if (!validation.passed) continue;
        
        cand.geometry = geom;
        cand.valid = false;
        return cand;
    }
    
    fmt::print("   ⚠️ Warning: could not generate valid candidate after {} attempts. Using base geometry.\n", MAX_ATTEMPTS);
    cand.geometry = base;
    return cand;
}

// ============================================================
// ВЫЧИСЛЕНИЕ ОШИБКИ (ТОЛЬКО СООТНОШЕНИЯ + ШТРАФЫ)
// ============================================================
double EvolutionaryOptimizer::computeError(const CCXData& ccx_data, const BellGeometry& geometry) {
    if (ccx_data.modes.size() < 5) {
        return 1e9;
    }
    
    const auto& modes = ccx_data.modes;
    
    double hum_freq = modes[0].frequency_hz;
    double hum_amp = modes[0].amplitude;
    double prime_freq = modes[1].frequency_hz;
    
    if (hum_freq < 1.0 || prime_freq < 1.0) {
        return 1e9;
    }
    
    auto sigmoid = [](double x, double k) -> double {
        if (x < -100.0) return 0.0;
        if (x > 100.0) return 1.0;
        return 2.0 / (1.0 + std::exp(-k * x)) - 1.0;
    };
    
    double freq_error = 0.0;
    int mode_indices[4] = {0, 2, 3, 4};
    
    for (int i = 0; i < 4; ++i) {
        int idx = mode_indices[i];
        if (idx >= (int)modes.size()) {
            return 1e9;
        }
        double freq = modes[idx].frequency_hz;
        double ratio = freq / prime_freq;
        double error_raw = targets_.freq_weights[i] * std::pow(ratio - targets_.freq_ratios[i], 2.0);
        freq_error += sigmoid(error_raw, targets_.k_freq);
    }
    
    double amp_error = 0.0;
    int amp_mode_indices[4] = {1, 2, 3, 4};
    
    for (int i = 0; i < 4; ++i) {
        int idx = amp_mode_indices[i];
        if (idx >= (int)modes.size()) {
            return 1e9;
        }
        double amp = modes[idx].amplitude;
        double amp_ratio = amp / hum_amp;
        double error_raw = targets_.amp_weights[i] * std::pow(amp_ratio - targets_.amp_ratios[i], 2.0);
        amp_error += sigmoid(error_raw, targets_.k_amp);
    }
    
    double angle_penalty = 0.0;
    const double max_angle_deg = 45.0;
    const double max_slope = std::tan(max_angle_deg * M_PI / 180.0);
    
    for (size_t i = 1; i < geometry.r_knots.size(); ++i) {
        double dz = geometry.getKnotZ(i) - geometry.getKnotZ(i-1);
        double dr = (geometry.r_knots[i] - geometry.r_knots[i-1]) * geometry.R_mouth;
        if (dz < 1e-9) continue;
        double slope = std::abs(dr / dz);
        if (slope > max_slope) {
            double excess = slope - max_slope;
            angle_penalty += excess * excess;
        }
    }
    angle_penalty /= (geometry.r_knots.size() - 1);
    double angle_error = 0.5 * sigmoid(angle_penalty, 2.0);
    
    double proportion_penalty = 0.0;
    double diameter = 2.0 * geometry.R_mouth;
    double ratio_diam = diameter / geometry.H;
    
    double target_min = 0.7;
    double target_max = 1.3;
    
    if (ratio_diam < target_min) {
        double excess = (target_min - ratio_diam) / target_min;
        proportion_penalty += excess * excess;
    }
    if (ratio_diam > target_max) {
        double excess = (ratio_diam - target_max) / target_max;
        proportion_penalty += excess * excess;
    }
    
    double proportion_error = 0.5 * sigmoid(proportion_penalty, 2.0);
    
    double total_error = freq_error + amp_error + angle_error + proportion_error;
    
    return total_error;
}

// ============================================================
// ОЦЕНКА КАНДИДАТА (С ПРОПОРЦИОНАЛЬНЫМИ ПРОВЕРКАМИ)
// ============================================================
double EvolutionaryOptimizer::evaluateCandidate(BellGeometry& geometry, CCXData& out_data) {
    BellMeshConfig config;
    config.nNodes = 15;
    
    BellGenerator generator(material_, config);
    generator.setGeometry(geometry);
    
    auto validation = generator.validateGeometry();
    if (!validation.passed) {
        return 1e9;
    }
    
    if (!generator.buildMesh()) {
        return 1e9;
    }
    
    if (!generator.saveToFile("bell.inp")) {
        return 1e9;
    }
    
#ifdef _WIN32
    int ccx_result = std::system("ccx bell > nul 2>&1");
#else
    int ccx_result = std::system("ccx bell > /dev/null 2>&1");
#endif
    if (ccx_result != 0) {
        return 1e9;
    }
    
    CCXData first_data = BellOptimizer::CCXParser::parse_file("bell.dat", true);
    if (first_data.modes.empty()) {
        return 1e9;
    }
    
    double prime_freq = 0.0;
    if (first_data.modes.size() > 1) {
        prime_freq = first_data.modes[1].frequency_hz;
    }
    
    if (prime_freq < 1.0) {
        return 1e9;
    }
    
    double scale = prime_freq / targets_.prime_freq_target;
    
    if (scale < 0.05 || scale > 20.0) {
        return 1e9;
    }
    
    BellGeometry scaled_geometry = geometry;
    scaled_geometry.scale(scale);
    
    // ============================================================
    // ПРОВЕРКИ ПОСЛЕ МАСШТАБИРОВАНИЯ (ТОЛЬКО ПРОПОРЦИИ)
    // ============================================================
    
    double t_top_ratio = scaled_geometry.t_top / scaled_geometry.H;
    double t_bottom_ratio = scaled_geometry.t_bottom / scaled_geometry.H;
    
    if (t_top_ratio < 0.005 || t_top_ratio > 0.050) {
        return 1e9;
    }
    if (t_bottom_ratio < 0.010 || t_bottom_ratio > 0.100) {
        return 1e9;
    }
    
    double diameter = 2.0 * scaled_geometry.R_mouth;
    double ratio_diam = diameter / scaled_geometry.H;
    if (ratio_diam < 0.7 || ratio_diam > 1.3) {
        return 1e9;
    }
    
    // ============================================================
    // ВТОРОЙ РАСЧЁТ
    // ============================================================
    generator.setGeometry(scaled_geometry);
    
    if (!generator.buildMesh()) {
        return 1e9;
    }
    
    if (!generator.saveToFile("bell.inp")) {
        return 1e9;
    }
    
#ifdef _WIN32
    ccx_result = std::system("ccx bell > nul 2>&1");
#else
    ccx_result = std::system("ccx bell > /dev/null 2>&1");
#endif
    if (ccx_result != 0) {
        return 1e9;
    }
    
    out_data = BellOptimizer::CCXParser::parse_file("bell.dat", true);
    if (out_data.modes.empty()) {
        return 1e9;
    }
    
    geometry = scaled_geometry;
    
    return computeError(out_data, scaled_geometry);
}

// ============================================================
// ИНИЦИАЛИЗАЦИЯ ПОПУЛЯЦИИ (С СОХРАНЕНИЕМ ЭЛИТЫ)
// ============================================================
std::vector<EvolutionaryOptimizer::Candidate> EvolutionaryOptimizer::initializePopulation() {
    std::vector<Candidate> population;
    population.reserve(evo_params_.population_size);
    
    if (has_initial_guess_) {
        fmt::print("📂 Using loaded geometry as initial guess\n");
        fmt::print("   Error: {:.6f}\n", best_error_);
        
        double max_error = 1.0;
        double scale = std::min(best_error_ / max_error, 1.0);
        double radius = 0.03 + 0.17 * scale;
        fmt::print("   Adaptive search radius: {:.0f}%\n", radius * 100);
        
        Candidate initial_cand;
        initial_cand.geometry = initial_guess_;
        initial_cand.error = best_error_;
        initial_cand.ccx_data = best_ccx_data_;
        initial_cand.valid = true;
        population.push_back(initial_cand);
        
        for (int i = 1; i < evo_params_.population_size; ++i) {
            Candidate cand = createNearCandidate(initial_guess_, best_error_);
            population.push_back(cand);
        }
        
        return population;
    }
    
    fmt::print("🆕 No saved geometry. Using reasonable shapes + randoms.\n");
    auto reasonable = getReasonableGeometries();
    for (const auto& geom : reasonable) {
        Candidate cand;
        cand.geometry = geom;
        cand.error = 1e9;
        cand.valid = false;
        population.push_back(cand);
    }
    
    while ((int)population.size() < evo_params_.population_size) {
        population.push_back(createRandomCandidate());
    }
    
    return population;
}

// ============================================================
// ОЦЕНКА ПОПУЛЯЦИИ
// ============================================================
void EvolutionaryOptimizer::evaluatePopulation(std::vector<Candidate>& population) {
    for (auto& cand : population) {
        if (cand.valid) continue;
        
        CCXData ccx_data;
        double error = evaluateCandidate(cand.geometry, ccx_data);
        
        cand.error = error;
        cand.ccx_data = ccx_data;
        cand.valid = true;
    }
}

// ============================================================
// ОТБОР
// ============================================================
std::vector<EvolutionaryOptimizer::Candidate> EvolutionaryOptimizer::selectElite(const std::vector<Candidate>& population) {
    std::vector<Candidate> sorted = population;
    std::sort(sorted.begin(), sorted.end(), 
              [](const Candidate& a, const Candidate& b) {
                  return a.error < b.error;
              });
    
    int num_elite = std::min(evo_params_.num_elite, (int)sorted.size());
    std::vector<Candidate> elite(sorted.begin(), sorted.begin() + num_elite);
    
    if (!elite.empty() && elite.front().error < 1e8) {
        if (elite.front().error < best_error_) {
            fmt::print("   🆕 New best! Error: {:.6f} (was {:.6f})\n", 
                       elite.front().error, best_error_);
            best_error_ = elite.front().error;
            best_geometry_ = elite.front().geometry;
            best_ccx_data_ = elite.front().ccx_data;
            saveBestResult();
        } else {
            fmt::print("   ℹ️  No improvement. Best: {:.6f}, Current: {:.6f}\n", 
                       best_error_, elite.front().error);
        }
    }
    
    return elite;
}

// ============================================================
// СКРЕЩИВАНИЕ
// ============================================================
EvolutionaryOptimizer::Candidate EvolutionaryOptimizer::crossover(const Candidate& p1, const Candidate& p2) {
    Candidate child;
    child.geometry = p1.geometry;
    
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int point = dist(rng_) * 3;
    
    if (point == 0) {
        child.geometry.H = p2.geometry.H;
    } else if (point == 1) {
        child.geometry.R_mouth = p2.geometry.R_mouth;
    } else if (point == 2) {
        for (size_t i = 1; i < child.geometry.r_knots.size() - 1; ++i) {
            child.geometry.r_knots[i] = p2.geometry.r_knots[i];
        }
    }
    
    if (dist(rng_) > 0.5) {
        child.geometry.t_top = p2.geometry.t_top;
    }
    if (dist(rng_) > 0.5) {
        child.geometry.t_bottom = p2.geometry.t_bottom;
    }
    
    if (!child.geometry.isValid()) {
        return createRandomCandidate();
    }
    
    return child;
}

// ============================================================
// МУТАЦИЯ
// ============================================================
void EvolutionaryOptimizer::mutate(BellGeometry& geom) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::uniform_real_distribution<double> distMutate(-1.0, 1.0);
    
    if (dist(rng_) < evo_params_.mutation_probability) {
        double delta = distMutate(rng_) * 0.05 * 1500.0;
        geom.H = std::max(500.0, std::min(2000.0, geom.H + delta));
    }
    
    if (dist(rng_) < evo_params_.mutation_probability) {
        double r_min = 0.4 * geom.H;
        double r_max = 0.8 * geom.H;
        double delta = distMutate(rng_) * 0.05 * (r_max - r_min);
        geom.R_mouth = std::max(r_min, std::min(r_max, geom.R_mouth + delta));
    }
    
    if (dist(rng_) < evo_params_.mutation_probability) {
        double current_ratio = geom.t_top / geom.H;
        double delta = distMutate(rng_) * 0.005;
        current_ratio = std::max(0.005, std::min(0.03, current_ratio + delta));
        geom.t_top = current_ratio * geom.H;
    }
    
    if (dist(rng_) < evo_params_.mutation_probability) {
        double current_ratio = geom.t_bottom / geom.H;
        double delta = distMutate(rng_) * 0.01;
        current_ratio = std::max(0.02, std::min(0.09, current_ratio + delta));
        geom.t_bottom = current_ratio * geom.H;
    }
    
    if (dist(rng_) < evo_params_.mutation_probability) {
        for (size_t i = 1; i < geom.r_knots.size() - 1; ++i) {
            double r_norm = geom.r_knots[i];
            double delta = distMutate(rng_) * 0.05;
            r_norm = std::max(0.0, std::min(1.0, r_norm + delta));
            geom.r_knots[i] = r_norm;
        }
    }
}

// ============================================================
// ЗАПУСК ОПТИМИЗАЦИИ
// ============================================================
void EvolutionaryOptimizer::run() {
    fmt::print("\n{:=^60}\n", "");
    fmt::print("{:=^60}\n", "  EVOLUTIONARY OPTIMIZER (Orthodox Bell Tuning)  ");
    fmt::print("{:=^60}\n", "");
    fmt::print("  Target Prime: {:.1f} Hz\n", targets_.prime_freq_target);
    fmt::print("  Population: {}\n", evo_params_.population_size);
    fmt::print("  Generations: {}\n", evo_params_.generations);
    fmt::print("  Mutation: {:.0f}%\n", evo_params_.mutation_probability * 100);
    fmt::print("{:=^60}\n\n", "");
    
    auto start_total = std::chrono::high_resolution_clock::now();
    
    std::vector<Candidate> population = initializePopulation();
    evaluatePopulation(population);
    
    auto elite = selectElite(population);
    
    if (!elite.empty() && elite.front().error < 1e8) {
        best_error_ = elite.front().error;
        best_geometry_ = elite.front().geometry;
        best_ccx_data_ = elite.front().ccx_data;
        saveBestResult();
    } else {
        fmt::print("⚠️ No valid candidates\n");
        best_geometry_ = getDefaultGeometry();
        best_error_ = 1e9;
    }
    
    printHeader();
    printGeneration(0, elite.empty() ? Candidate() : elite.front());
    
    for (int gen = 1; gen <= evo_params_.generations; ++gen) {
        fmt::print("\n--- Generation {} ---\n", gen);
        
        std::vector<Candidate> new_population;
        new_population.reserve(evo_params_.population_size);
        
        for (const auto& c : elite) {
            new_population.push_back(c);
        }
        
        while ((int)new_population.size() < evo_params_.population_size) {
            std::uniform_int_distribution<int> elite_dist(0, elite.size() - 1);
            const auto& p1 = elite[elite_dist(rng_)];
            const auto& p2 = elite[elite_dist(rng_)];
            
            Candidate child = crossover(p1, p2);
            mutate(child.geometry);
            new_population.push_back(child);
        }
        
        evaluatePopulation(new_population);
        elite = selectElite(new_population);
        
        printGeneration(gen, elite.front());
        
        if (best_error_ < 0.01) {
            fmt::print("\n🎯 Target achieved!\n");
            break;
        }
        
        population = new_population;
    }
    
    fmt::print("\n📝 Final save...\n");
    if (!elite.empty() && elite.front().error < 1e8) {
        if (elite.front().error < best_error_) {
            best_error_ = elite.front().error;
            best_geometry_ = elite.front().geometry;
            best_ccx_data_ = elite.front().ccx_data;
        }
        saveBestResult();
    }
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count();
    
    fmt::print("\n{:=^60}\n", "");
    fmt::print("{:=^60}\n", "  OPTIMIZATION COMPLETED  ");
    fmt::print("{:=^60}\n", "");
    fmt::print("  Best error: {:.6f}\n", best_error_);
    fmt::print("  Time: {:.1f} s\n", duration / 1000.0);
    fmt::print("{:=^60}\n", "");
    
    printFinal(elite.front());
}

// ============================================================
// ВЫВОД РЕЗУЛЬТАТОВ
// ============================================================
void EvolutionaryOptimizer::printHeader() {
    fmt::print("\n");
    fmt::print("{:>4} {:>12} {:>10} {:>10} {:>10} {:>10} {:>12} {:>12} {:>12} {:>12} {:>12}\n",
               "Gen", "Error", "H", "Rm", "t_top", "t_bot", 
               "Hum", "Prime", "Tertia", "Quinta", "Nominal");
    fmt::print("{:-<120}\n", "");
}

void EvolutionaryOptimizer::printGeneration(int gen, const Candidate& best) {
    if (best.error >= 1e8) {
        fmt::print("{:>4} {:>12}\n", gen, "FAILED");
        return;
    }
    
    const auto& modes = best.ccx_data.modes;
    
    double hum_freq = 0.0, prime_freq = 0.0, tertia_freq = 0.0, quinta_freq = 0.0, nominal_freq = 0.0;
    
    if (modes.size() > 0) hum_freq = modes[0].frequency_hz;
    if (modes.size() > 1) prime_freq = modes[1].frequency_hz;
    if (modes.size() > 2) tertia_freq = modes[2].frequency_hz;
    if (modes.size() > 3) quinta_freq = modes[3].frequency_hz;
    if (modes.size() > 4) nominal_freq = modes[4].frequency_hz;
    
    fmt::print("{:>4} {:>12.4f} {:>10.1f} {:>10.1f} {:>10.1f} {:>10.1f} {:>12.1f} {:>12.1f} {:>12.1f} {:>12.1f} {:>12.1f}\n",
               gen, best.error, best.geometry.H, best.geometry.R_mouth,
               best.geometry.t_top, best.geometry.t_bottom,
               hum_freq, prime_freq, tertia_freq, quinta_freq, nominal_freq);
}

void EvolutionaryOptimizer::printFinal(const Candidate& best) {
    if (best.error >= 1e8) {
        fmt::print("\n  No valid solution found.\n");
        return;
    }
    
    fmt::print("\n  Best geometry:\n");
    fmt::print("    H         = {:.1f} mm\n", best.geometry.H);
    fmt::print("    R_mouth   = {:.1f} mm\n", best.geometry.R_mouth);
    fmt::print("    t_top     = {:.1f} mm\n", best.geometry.t_top);
    fmt::print("    t_bottom  = {:.1f} mm\n", best.geometry.t_bottom);
    fmt::print("    Radii     = ");
    for (double r : best.geometry.r_knots) {
        fmt::print("{:.3f} ", r);
    }
    fmt::print("\n");
    
    if (best.ccx_data.modes.size() >= 5) {
        fmt::print("\n  Frequencies:\n");
        const char* names[] = {"Hum", "Prime", "Tertia", "Quinta", "Nominal"};
        for (int i = 0; i < 5 && i < (int)best.ccx_data.modes.size(); ++i) {
            fmt::print("    {}: {:.2f} Hz", names[i], best.ccx_data.modes[i].frequency_hz);
            if (i > 0) {
                double ratio = best.ccx_data.modes[i].frequency_hz / best.ccx_data.modes[1].frequency_hz;
                fmt::print(" (ratio: {:.3f})", ratio);
            }
            fmt::print("\n");
        }
    }
}