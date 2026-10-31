// bell_optimizer.cpp
#include "bell_optimizer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "simple_json.hpp"
#include "../excitation/meridional_influence.hpp"
#include "../math/modal_coordinate_morphing.hpp"
#include "../fem/cylindrical_modal_evaluator.hpp"

// ============================================================
// СТРУКТУРА ДЛЯ СОХРАНЕНИЯ РЕЗУЛЬТАТА
// ============================================================
struct BestResult {
    double error = 0.0;
    BellGeometry geometry;
    ModalData modal_data;
    
    bool save(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return false;

        file << "{\n";
        file << "    \"error\": " << error << ",\n";
        file << "    \"frequencies\": {\n";
        if (modal_data.modes.size() > 0) file << "        \"hum\": " << modal_data.modes[0].frequency_hz << ",\n";
        if (modal_data.modes.size() > 1) file << "        \"prime\": " << modal_data.modes[1].frequency_hz << ",\n";
        if (modal_data.modes.size() > 2) file << "        \"tertia\": " << modal_data.modes[2].frequency_hz << ",\n";
        if (modal_data.modes.size() > 3) file << "        \"quinta\": " << modal_data.modes[3].frequency_hz << ",\n";
        if (modal_data.modes.size() > 4) file << "        \"nominal\": " << modal_data.modes[4].frequency_hz << "\n";
        file << "    },\n";
        file << "    \"geometry\": {\n";
        file << "        \"H\": " << geometry.H << ",\n";
        file << "        \"R_mouth\": " << geometry.R_mouth << ",\n";
        file << "        \"r_knots\": [";
        for (std::size_t i = 0; i < geometry.r_knots.size(); ++i) {
            if (i != 0) file << ", ";
            file << geometry.r_knots[i];
        }
        file << "],\n";
        file << "        \"t_bottom\": " << geometry.t_bottom << ",\n";
        file << "        \"t_top\": " << geometry.t_top << ",\n";
        file << "        \"ellipticity\": " << geometry.ellipticity << "\n";
        file << "    }\n";
        file << "}\n";
        return file.good();
    }

    bool load(const std::string& filename) {
        try {
            const auto data = bell::json::parse_file(filename);
            error = bell::json::number(data, "error", 0.0);
            const auto* g = data.find("geometry");
            if (g && g->is_object()) {
                geometry.H = bell::json::number(*g, "H", 1000.0);
                geometry.R_mouth = bell::json::number(*g, "R_mouth", 400.0);
                geometry.t_top = bell::json::number(*g, "t_top", 10.0);
                geometry.t_bottom = bell::json::number(*g, "t_bottom", 60.0);
                geometry.ellipticity = bell::json::number(*g, "ellipticity", 1.0);
                if (const auto* knots = g->find("r_knots"); knots && knots->is_array()) {
                    geometry.r_knots.clear();
                    for (const auto& item : knots->array()) {
                        geometry.r_knots.push_back(item.number(0.0));
                    }
                }
            }
            return true;
        } catch (const std::exception& e) {
            std::cerr << "❌ Ошибка загрузки JSON: " << e.what() << "\n";
            return false;
        }
    }

};

// ============================================================
// КОНСТРУКТОР
// ============================================================
EvolutionaryOptimizer::EvolutionaryOptimizer(const BellMaterial& material,
                                             const TargetData& targets,
                                             const EvolutionParams& evo_params,
                                             std::unique_ptr<ModalSolver> solver,
                                             std::unique_ptr<Models::GeometrySeedProvider> seed_provider)
    : material_(material)
    , targets_(targets)
    , evo_params_(evo_params)
    , seed_provider_(seed_provider ? std::move(seed_provider) : Models::make_seed_provider("neutral"))
    , rng_(evo_params_.deterministic ? static_cast<std::mt19937::result_type>(evo_params_.seed) : rd_())
    , solver_(std::move(solver))
{
    if (!solver_) throw std::invalid_argument("EvolutionaryOptimizer requires a ModalSolver");
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
    if (!seed_provider_) return {};
    return seed_provider_->make_seeds(material_, evo_params_.num_knots);
}

// ============================================================
// УСТАНОВИТЬ НАЧАЛЬНОЕ ПРИБЛИЖЕНИЕ
// ============================================================
void EvolutionaryOptimizer::setInitialGuess(const BellGeometry& initial) {
    initial_guess_ = initial;
    if (!evo_params_.optimize_ellipticity) {
        initial_guess_.ellipticity = 1.0;
    } else {
        initial_guess_.ellipticity = std::clamp(initial_guess_.ellipticity,
            evo_params_.ellipticity_min, evo_params_.ellipticity_max);
    }
    has_initial_guess_ = true;
}

// The modal ratios and normalized participation ratios are invariant under
// uniform geometric scaling. Searching H as an independent variable therefore
// adds a redundant degree of freedom. The optimizer searches a dimensionless
// shape and restores the physical scale from the target Prime frequency.
BellGeometry EvolutionaryOptimizer::canonicalizeShape(const BellGeometry& geometry, double reference_height) {
    BellGeometry shape = geometry;
    if (!(reference_height > 0.0) || !(geometry.H > 0.0) || !std::isfinite(geometry.H)) {
        return shape;
    }
    const double scale = reference_height / geometry.H;
    shape.scale(scale);
    return shape;
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
    if (!evo_params_.optimize_ellipticity) {
        loaded_geometry.ellipticity = 1.0;
    } else {
        loaded_geometry.ellipticity = std::clamp(loaded_geometry.ellipticity,
            evo_params_.ellipticity_min, evo_params_.ellipticity_max);
    }
    
    std::cout << "📂 Loaded geometry from " << filename << "\n";
    std::cout << std::fixed << std::setprecision(1) << "   Old H = " << loaded_geometry.H << ", R_mouth = " << loaded_geometry.R_mouth << "\n";
    std::cout << std::fixed << std::setprecision(6) << "   Old error = " << result.error << "\n";
    std::cout << std::fixed << std::setprecision(1) << "   🔄 Rescaling to target frequency: " << targets_.prime_freq_target << " Hz\n";
    
    // ============================================================
    // 1. РАСЧИТЫВАЕМ ЗАГРУЖЕННУЮ ГЕОМЕТРИЮ (ОДИН РАЗ)
    // ============================================================
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material_, config);
    generator.setGeometry(loaded_geometry);
    
    if (!generator.buildMesh()) {
        std::cout << "   ❌ Failed to build mesh\n";
        return false;
    }
    
    Excitation::MeridionalInfluence influence;
    ModalProblem problem{generator.getMeshData(), material_, ".", "bell", &influence};
    ModalData first_data;
    if (!solver_->solve(problem, first_data)) {
        std::cout << "   ❌ Modal solver failed (" << solver_->name() << ")\n";
        return false;
    }
    if (evo_params_.optimize_ellipticity && loaded_geometry.ellipticity < 1.0) {
        const Math::CoordinateScale3D map{1.0, loaded_geometry.ellipticity, 1.0};
        for (auto& mode : first_data.modes) {
            if (mode.nodal_radial.empty() || mode.nodal_axial.empty()) {
                std::cout << "   ❌ Ellipticity search needs captured modal fields\n";
                return false;
            }
            const auto morph = Math::evaluate_mechanical_morphing(
                generator.getMeshData(), material_, mode, map);
            mode.frequency_hz = morph.morphed_frequency_hz;
            mode.eigenvalue = std::pow(2.0 * 3.14159265358979323846 * mode.frequency_hz, 2.0);
        }
        first_data.sort_by_frequency();
    }
    
    double prime_freq = 0.0;
    if (first_data.modes.size() > 1) {
        prime_freq = first_data.modes[1].frequency_hz;
    }
    
    if (prime_freq < 1.0) {
        std::cout << std::fixed << std::setprecision(2) << "   ❌ Invalid prime frequency: " << prime_freq << "\n";
        return false;
    }
    
    std::cout << std::fixed << std::setprecision(2) << "   Current Prime: " << prime_freq << " Hz\n";
    
    // ============================================================
    // 2. ВЫЧИСЛЯЕМ МАСШТАБ (ПРАВИЛЬНАЯ ФОРМУЛА!)
    // ============================================================
    double scale = prime_freq / targets_.prime_freq_target;
    
    if (scale < 0.05 || scale > 20.0) {
        std::cout << std::fixed << std::setprecision(3) << "   ❌ Invalid scale: " << scale << "\n";
        return false;
    }
    
    std::cout << std::fixed << std::setprecision(3) << "   Scale factor: " << scale << "\n";
    
    // ============================================================
    // 3. МАСШТАБИРУЕМ ГЕОМЕТРИЮ
    // ============================================================
    BellGeometry scaled_geometry = loaded_geometry;
    scaled_geometry.scale(scale);
    
    // ============================================================
    // 4. МАСШТАБИРУЕМ ЧАСТОТЫ АНАЛИТИЧЕСКИ (БЕЗ ВТОРОГО ЗАПУСКА!)
    // ============================================================
    ModalData scaled_data = first_data;
    for (auto& mode : scaled_data.modes) {
        mode.frequency_hz /= scale;  // ω_new = ω_old / scale
    }
    
    // ============================================================
    // 5. ПРОВЕРКИ ПОСЛЕ МАСШТАБИРОВАНИЯ (ПРОПОРЦИИ НЕ МЕНЯЮТСЯ!)
    // ============================================================
    double t_top_ratio = scaled_geometry.t_top / scaled_geometry.H;
    double t_bottom_ratio = scaled_geometry.t_bottom / scaled_geometry.H;
    
    if (t_top_ratio < 0.005 || t_top_ratio > 0.050) {
        std::cout << std::fixed << std::setprecision(4) << "   ❌ t_top/H ratio out of range: " << t_top_ratio << "\n";
        return false;
    }
    if (t_bottom_ratio < 0.010 || t_bottom_ratio > 0.100) {
        std::cout << std::fixed << std::setprecision(4) << "   ❌ t_bottom/H ratio out of range: " << t_bottom_ratio << "\n";
        return false;
    }
    
    double diameter = 2.0 * scaled_geometry.R_mouth;
    double ratio_diam = diameter / scaled_geometry.H;
    if (ratio_diam < 0.7 || ratio_diam > 1.3) {
        std::cout << std::fixed << std::setprecision(3) << "   ❌ Diameter/Height ratio out of range: " << ratio_diam << "\n";
        return false;
    }
    
    std::cout << std::fixed << std::setprecision(1) << "   Scaled H = " << scaled_geometry.H << ", R_mouth = " << scaled_geometry.R_mouth << "\n";
    std::cout << std::fixed << std::setprecision(4) << "   t_top/H = " << t_top_ratio << ", t_bottom/H = " << t_bottom_ratio << "\n";
    
    // ============================================================
    // 6. СОХРАНЯЕМ ПЕРЕСЧИТАННЫЕ ДАННЫЕ
    // ============================================================
    setInitialGuess(scaled_geometry);
    best_error_ = computeError(scaled_data, scaled_geometry);
    best_geometry_ = scaled_geometry;
    best_modal_data_ = scaled_data;
    modal_reference_ = scaled_data;
    has_modal_reference_ = !modal_reference_.modes.empty();
    
    std::cout << "✅ Rescaled successfully!\n";
    std::cout << std::fixed << std::setprecision(6) << "   New error: " << best_error_ << "\n";
    std::cout << std::fixed << std::setprecision(2) << "   New Prime: " << best_modal_data_.modes[1].frequency_hz << " Hz\n";
    
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
    std::cout << "   🔍 Saving best result...\n";
    std::cout << std::fixed << std::setprecision(6) << "      Error: " << best_error_ << "\n";
    std::cout << "      Modes: " << best_modal_data_.modes.size() << "\n";
    
    if (best_error_ >= 1e8) {
        std::cout << "      ❌ Error too high (>= 1e8)\n";
        return;
    }
    
    if (best_modal_data_.modes.size() < 3) {
        std::cout << "      ❌ Not enough modes: " << best_modal_data_.modes.size() << " < 3\n";
        return;
    }
    
    bool has_valid_freq = false;
    for (const auto& mode : best_modal_data_.modes) {
        if (mode.frequency_hz > 20.0) {
            has_valid_freq = true;
            break;
        }
    }
    
    if (!has_valid_freq) {
        std::cout << "      ❌ No valid frequencies (> 20 Hz)\n";
        return;
    }
    
    BestResult result;
    result.error = best_error_;
    result.geometry = best_geometry_;
    result.modal_data = best_modal_data_;
    
    std::cout << "      💾 Writing to best_geometry.json...\n";
    if (result.save("best_geometry.json")) {
        std::cout << "      ✅ Saved successfully!\n";
    } else {
        std::cout << "      ❌ SAVE FAILED! Cannot open file for writing.\n";
        std::cout << "      ⚠️  Check file permissions or if file is open in another program.\n";
    }
    
    if (!best_modal_data_.source_file.empty()) {
        copyFile(best_modal_data_.source_file, "best_bell.dat");
        std::cout << "      ✅ Copied modal result -> best_bell.dat\n";
    }
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
    std::uniform_real_distribution<double> distR(0.0, 1.0);
    
    BellGeometry geom;
    if (evo_params_.optimize_ellipticity) {
        std::uniform_real_distribution<double> qdist(evo_params_.ellipticity_min, evo_params_.ellipticity_max);
        geom.ellipticity = qdist(rng_);
    } else {
        geom.ellipticity = 1.0;
    }
    // Shape search is dimensionless; absolute size is recovered after the
    // modal calculation from the requested Prime frequency.
    geom.H = 1000.0;
    
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
        if (evo_params_.optimize_ellipticity) {
            geom.ellipticity = std::clamp(base.ellipticity + dist(rng_) * 0.10,
                                         evo_params_.ellipticity_min, evo_params_.ellipticity_max);
        } else {
            geom.ellipticity = 1.0;
        }
        
        // H is deliberately not mutated: it is a pure scale degree of freedom.
        geom.H = 1000.0;
        const double base_scale = (base.H > 0.0) ? (1000.0 / base.H) : 1.0;
        geom.R_mouth = base.R_mouth * base_scale;
        geom.t_top = base.t_top * base_scale;
        geom.t_bottom = base.t_bottom * base_scale;
        double r_min = 0.4 * geom.H;
        double r_max = 0.8 * geom.H;
        double delta_R = dist(rng_) * geom.R_mouth;
        geom.R_mouth = std::max(r_min, std::min(r_max, base.R_mouth + delta_R));
        
        double delta_t_top = dist(rng_) * geom.t_top;
        geom.t_top = std::max(2.0, std::min(60.0, base.t_top + delta_t_top));
        
        double delta_t_bottom = dist(rng_) * geom.t_bottom;
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
    
    std::cout << "   ⚠️ Warning: could not generate valid candidate after " << MAX_ATTEMPTS << " attempts. Using base geometry.\n";
    cand.geometry = base;
    return cand;
}

// ============================================================
// ВЫЧИСЛЕНИЕ ОШИБКИ (ТОЛЬКО СООТНОШЕНИЯ + ШТРАФЫ)
// ============================================================
double EvolutionaryOptimizer::computeError(const ModalData& modal_data, const BellGeometry& geometry) {
    if (modal_data.modes.size() < 5) return 1e9;

    const auto& modes = modal_data.modes;
    const double hum_freq = modes[0].frequency_hz;
    const double hum_amp = modes[0].amplitude;
    const double prime_freq = modes[1].frequency_hz;
    if (hum_freq < BellConstants::MIN_FREQ_HZ || prime_freq < BellConstants::MIN_FREQ_HZ ||
        hum_amp <= 1e-12) return 1e9;

    // Modal ratios are positive quantities. Log-relative error treats a +10%
    // and -10% deviation symmetrically and avoids the saturation of the old
    // sigmoid penalty.
    double freq_error = 0.0;
    const int freq_indices[4] = {0, 2, 3, 4};
    for (int i = 0; i < 4; ++i) {
        const int idx = freq_indices[i];
        if (idx >= static_cast<int>(modes.size()) || targets_.freq_ratios[i] <= 0.0) return 1e9;
        const double ratio = modes[idx].frequency_hz / prime_freq;
        if (!(ratio > 0.0) || !std::isfinite(ratio)) return 1e9;
        const double e = std::log(ratio / targets_.freq_ratios[i]);
        freq_error += targets_.freq_weights[i] * e * e;
    }

    double amp_error = 0.0;
    const int amp_indices[4] = {1, 2, 3, 4};
    for (int i = 0; i < 4; ++i) {
        const int idx = amp_indices[i];
        if (idx >= static_cast<int>(modes.size()) || targets_.amp_ratios[i] <= 0.0) return 1e9;
        const double ratio = modes[idx].amplitude / hum_amp;
        if (!(ratio > 0.0) || !std::isfinite(ratio)) return 1e9;
        const double e = std::log(ratio / targets_.amp_ratios[i]);
        amp_error += targets_.amp_weights[i] * e * e;
    }

    // The geometry generator performs the authoritative derivative-based slope
    // check. Keep the objective only responsible for smooth proportion penalty.
    double proportion_error = 0.0;
    const double diameter_ratio = 2.0 * geometry.R_mouth / geometry.H;
    if (diameter_ratio < BellConstants::DIAMETER_RATIO_MIN) {
        const double e = std::log(BellConstants::DIAMETER_RATIO_MIN / diameter_ratio);
        proportion_error += e * e;
    } else if (diameter_ratio > BellConstants::DIAMETER_RATIO_MAX) {
        const double e = std::log(diameter_ratio / BellConstants::DIAMETER_RATIO_MAX);
        proportion_error += e * e;
    }

    return freq_error + amp_error + proportion_error;
}

// ============================================================
// ОЦЕНКА КАНДИДАТА (ОДИН ЗАПУСК + АНАЛИТИЧЕСКОЕ МАСШТАБИРОВАНИЕ)
// ============================================================
double EvolutionaryOptimizer::evaluateCandidate(BellGeometry& geometry, ModalData& out_data) {
    geometry = canonicalizeShape(geometry);
    BellMeshConfig config;
    config.nNodes = 15;
    
    // ============================================================
    // ЭТАП 1: ГЕНЕРАЦИЯ СЕТКИ
    // ============================================================
    BellGenerator generator(material_, config);
    generator.setGeometry(geometry);
    
    auto validation = generator.validateGeometry();
    if (!validation.passed) {
        return 1e9;
    }
    
    if (!generator.buildMesh()) {
        return 1e9;
    }
    
    Excitation::MeridionalInfluence influence;
    ModalProblem problem{generator.getMeshData(), material_, ".", "bell", &influence};
    
    // ============================================================
    // ЭТАП 2: ОДИН ЗАПУСК МОДАЛЬНОГО РЕШАТЕЛЯ
    // ============================================================
    ModalData raw_data;
    if (!solver_->solve(problem, raw_data)) {
        return 1e9;
    }

    // Modal identity must not be inferred from frequency rank alone. During
    // geometry optimization two eigenfrequencies can cross; the physical
    // family is then better represented by the displacement shape. The
    // reference is established from the first valid solve (or from a loaded
    // best result) and remains fixed during the run.
    if (!has_modal_reference_) {
        modal_reference_ = raw_data;
        has_modal_reference_ = !modal_reference_.modes.empty();
    } else {
        const auto tracking = trackModalFamilies(modal_reference_, raw_data);
        if (tracking.used_shape_tracking && tracking.ambiguous) {
            // Do not reject a physically valid candidate merely because a pair
            // is near-degenerate. Record the diagnostic and let the objective
            // continue with the best available assignment.
            std::cout << "   ⚠️ Modal-family tracking: near-degenerate/ambiguous pair, "
                      << "min shape similarity=" << std::fixed << std::setprecision(3)
                      << tracking.minimum_similarity << "\n" << std::defaultfloat;
        }
    }
    
    // ============================================================
    // ЭТАП 3: РЕЗУЛЬТАТ МОДАЛЬНОГО РЕШЕЛЯ
    // ============================================================
    if (raw_data.modes.empty()) {
        return 1e9;
    }

    // Ellipticity is an optional search degree of freedom. Circular-bell runs
    // pin q=1 exactly. For an enabled elliptical search the direct vector FEM
    // path can replace the historical Rayleigh coordinate morph. The legacy
    // morph remains available as an explicit fallback until the Chinese-bell
    // objective is finalized and promoted to production default.
    if (!evo_params_.optimize_ellipticity) {
        geometry.ellipticity = 1.0;
    } else if (geometry.ellipticity < evo_params_.ellipticity_min ||
               geometry.ellipticity > evo_params_.ellipticity_max) {
        return 1e9;
    } else if (geometry.ellipticity < 1.0 && evo_params_.use_direct_elliptical_fem) {
        try {
            const auto fem_model = FEM::make_cylindrical_bell_model(geometry, material_, config);
            FEM::CylindricalModalEvaluationConfig fem_config;
            fem_config.harmonics = {0, 2, 4};
            fem_config.angular_samples = 64;
            fem_config.requested_modes = 10;
            fem_config.max_iterations = 1200;
            fem_config.tolerance = 1.0e-6;
            fem_config.minimum_frequency_hz = BellConstants::MIN_FREQ_HZ;
            Excitation::MeridionalInfluence fem_influence;
            ModalData direct_data;
            if (!FEM::evaluate_cylindrical_bell_modes(
                    fem_model, fem_influence.make(fem_model.mesh_si), fem_config, direct_data)) {
                return 1e9;
            }
            raw_data = std::move(direct_data);
        } catch (const std::exception&) {
            return 1e9;
        }
    } else if (geometry.ellipticity < 1.0) {
        const Math::CoordinateScale3D map{1.0, geometry.ellipticity, 1.0};
        for (auto& mode : raw_data.modes) {
            if (mode.nodal_radial.empty() || mode.nodal_axial.empty()) return 1e9;
            const auto morph = Math::evaluate_mechanical_morphing(
                generator.getMeshData(), material_, mode, map);
            mode.frequency_hz = morph.morphed_frequency_hz;
            mode.eigenvalue = std::pow(2.0 * 3.14159265358979323846 * mode.frequency_hz, 2.0);
        }
        raw_data.sort_by_frequency();
    }
    
    double prime_freq = 0.0;
    if (raw_data.modes.size() > 1) {
        prime_freq = raw_data.modes[1].frequency_hz;
    }
    
    if (prime_freq < 1.0) {
        return 1e9;
    }
    
    // ============================================================
    // ЭТАП 4: ВЫЧИСЛЯЕМ МАСШТАБ
    // ============================================================
    double scale = prime_freq / targets_.prime_freq_target;
    
    if (scale < 0.05 || scale > 20.0) {
        return 1e9;
    }
    
    // ============================================================
    // ЭТАП 5: МАСШТАБИРУЕМ ГЕОМЕТРИЮ (ДЛЯ СОХРАНЕНИЯ)
    // ============================================================
    geometry.scale(scale);
    
    // ============================================================
    // ЭТАП 6: МАСШТАБИРУЕМ ЧАСТОТЫ АНАЛИТИЧЕСКИ
    // ============================================================
    out_data = raw_data;
    for (auto& mode : out_data.modes) {
        mode.frequency_hz /= scale;  // ω_new = ω_old / scale
    }
    
    // ============================================================
    // ЭТАП 7: ПРОВЕРКИ ПРОПОРЦИЙ (НЕ МЕНЯЮТСЯ ПРИ МАСШТАБЕ!)
    // ============================================================
    double t_top_ratio = geometry.t_top / geometry.H;
    double t_bottom_ratio = geometry.t_bottom / geometry.H;
    
    if (t_top_ratio < 0.005 || t_top_ratio > 0.050) {
        return 1e9;
    }
    if (t_bottom_ratio < 0.010 || t_bottom_ratio > 0.100) {
        return 1e9;
    }
    
    double diameter = 2.0 * geometry.R_mouth;
    double ratio_diam = diameter / geometry.H;
    if (ratio_diam < 0.7 || ratio_diam > 1.3) {
        return 1e9;
    }
    
    // ============================================================
    // ЭТАП 8: ВЫЧИСЛЯЕМ ОШИБКУ
    // ============================================================
    return computeError(out_data, geometry);
}

// ============================================================
// ИНИЦИАЛИЗАЦИЯ ПОПУЛЯЦИИ (С СОХРАНЕНИЕМ ЭЛИТЫ)
// ============================================================
std::vector<EvolutionaryOptimizer::Candidate> EvolutionaryOptimizer::initializePopulation() {
    std::vector<Candidate> population;
    population.reserve(evo_params_.population_size);
    
    if (has_initial_guess_) {
        std::cout << "📂 Using loaded geometry as initial guess\n";
        std::cout << std::fixed << std::setprecision(6) << "   Error: " << best_error_ << "\n";
        
        double max_error = 1.0;
        double scale = std::min(best_error_ / max_error, 1.0);
        double radius = 0.03 + 0.17 * scale;
        std::cout << std::fixed << std::setprecision(0) << "   Adaptive search radius: " << radius * 100 << "%\n";
        
        Candidate initial_cand;
        initial_cand.geometry = initial_guess_;
        initial_cand.error = best_error_;
        initial_cand.modal_data = best_modal_data_;
        initial_cand.valid = true;
        population.push_back(initial_cand);
        
        for (int i = 1; i < evo_params_.population_size; ++i) {
            Candidate cand = createNearCandidate(initial_guess_, best_error_);
            population.push_back(cand);
        }
        
        return population;
    }
    
    std::cout << "🆕 No saved geometry. Using reasonable shapes + randoms.\n";
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
        
        ModalData modal_data;
        double error = evaluateCandidate(cand.geometry, modal_data);
        
        cand.error = error;
        cand.modal_data = modal_data;
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
            std::cout << std::fixed << std::setprecision(6) << "   🆕 New best! Error: " << elite.front().error << " (was " << best_error_ << ")\n";
            best_error_ = elite.front().error;
            best_geometry_ = elite.front().geometry;
            best_modal_data_ = elite.front().modal_data;
            saveBestResult();
        } else {
            std::cout << std::fixed << std::setprecision(6) << "   ℹ️  No improvement. Best: " << best_error_ << ", Current: " << elite.front().error << "\n";
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
    int point = dist(rng_) * (evo_params_.optimize_ellipticity ? 4 : 3);
    
    if (point == 0) {
        child.geometry.H = p2.geometry.H;
    } else if (point == 1) {
        child.geometry.R_mouth = p2.geometry.R_mouth;
    } else if (point == 2) {
        for (size_t i = 1; i < child.geometry.r_knots.size() - 1; ++i) {
            child.geometry.r_knots[i] = p2.geometry.r_knots[i];
        }
    } else if (point == 3 && evo_params_.optimize_ellipticity) {
        child.geometry.ellipticity = p2.geometry.ellipticity;
    }
    if (!evo_params_.optimize_ellipticity) child.geometry.ellipticity = 1.0;
    
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
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    const double s = std::clamp(evo_params_.mutation_strength, 0.0, 1.0);

    // H is a scale variable and is intentionally excluded from the shape
    // search. The following perturbations operate on dimensionless geometry.
    // Circular profiles freeze q exactly; ellipticity is searched only by explicit opt-in.
    if (evo_params_.optimize_ellipticity && dist(rng_) < evo_params_.mutation_probability) {
        geom.ellipticity = std::clamp(geom.ellipticity + unit(rng_) * s * 0.15,
                                      evo_params_.ellipticity_min, evo_params_.ellipticity_max);
    } else if (!evo_params_.optimize_ellipticity) {
        geom.ellipticity = 1.0;
    }
    if (dist(rng_) < evo_params_.mutation_probability) {
        const double ratio = std::clamp(geom.R_mouth / geom.H + unit(rng_) * s * 0.10,
                                        BellConstants::MIN_RMOUTH_RATIO,
                                        BellConstants::MAX_RMOUTH_RATIO);
        geom.R_mouth = ratio * geom.H;
    }

    if (dist(rng_) < evo_params_.mutation_probability) {
        const double ratio = std::clamp(geom.t_top / geom.H + unit(rng_) * s * 0.01,
                                        BellConstants::MIN_TTOP_RATIO,
                                        BellConstants::MAX_TTOP_RATIO);
        geom.t_top = ratio * geom.H;
    }

    if (dist(rng_) < evo_params_.mutation_probability) {
        const double ratio = std::clamp(geom.t_bottom / geom.H + unit(rng_) * s * 0.02,
                                        BellConstants::MIN_TBOTTOM_RATIO,
                                        BellConstants::MAX_TBOTTOM_RATIO);
        geom.t_bottom = ratio * geom.H;
    }

    if (dist(rng_) < evo_params_.mutation_probability) {
        for (size_t i = 1; i + 1 < geom.r_knots.size(); ++i) {
            geom.r_knots[i] = std::clamp(geom.r_knots[i] + unit(rng_) * s * 0.10,
                                          0.0, 1.0);
        }
        // Projection onto the feasible monotone knot set.
        for (size_t i = 1; i + 1 < geom.r_knots.size(); ++i) {
            geom.r_knots[i] = std::max(geom.r_knots[i], geom.r_knots[i - 1] + 1e-4);
        }
        for (size_t i = geom.r_knots.size() - 1; i > 1; --i) {
            geom.r_knots[i - 1] = std::min(geom.r_knots[i - 1], geom.r_knots[i] - 1e-4);
        }
    }
}


// ============================================================
// ЗАПУСК ОПТИМИЗАЦИИ
// ============================================================
void EvolutionaryOptimizer::run() {
    std::cout << "\n============================================================\n  EVOLUTIONARY OPTIMIZER (Orthodox Bell Tuning)  \n";
    std::cout << "============================================================\n";
    std::cout << std::fixed << std::setprecision(1) << "  Target Prime: " << targets_.prime_freq_target << " Hz\n";
    std::cout << "  Population: " << evo_params_.population_size << "\n";
    std::cout << "  Generations: " << evo_params_.generations << "\n";
    std::cout << std::fixed << std::setprecision(0) << "  Mutation: " << evo_params_.mutation_probability * 100 << "%\n";
    std::cout << "  Mode: one run per candidate via pluggable modal solver (analytical scaling)\n";
    std::cout << "  Ellipticity search: " << (evo_params_.optimize_ellipticity ? "enabled" : "frozen at q=1") << "\n";
    std::cout << "============================================================\n\n";
    
    auto start_total = std::chrono::high_resolution_clock::now();
    
    std::vector<Candidate> population = initializePopulation();
    evaluatePopulation(population);
    
    auto elite = selectElite(population);
    
    if (!elite.empty() && elite.front().error < 1e8) {
        best_error_ = elite.front().error;
        best_geometry_ = elite.front().geometry;
        best_modal_data_ = elite.front().modal_data;
        saveBestResult();
    } else {
        std::cout << "⚠️ No valid candidates\n";
        best_geometry_ = getDefaultGeometry();
        best_error_ = 1e9;
    }
    
    printHeader();
    printGeneration(0, elite.empty() ? Candidate() : elite.front());
    
    for (int gen = 1; gen <= evo_params_.generations; ++gen) {
        std::cout << "\n--- Generation " << gen << " ---\n";
        
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
            std::cout << "\n🎯 Target achieved!\n";
            break;
        }
        
        population = new_population;
    }
    
    std::cout << "\n📝 Final save...\n";
    if (!elite.empty() && elite.front().error < 1e8) {
        if (elite.front().error < best_error_) {
            best_error_ = elite.front().error;
            best_geometry_ = elite.front().geometry;
            best_modal_data_ = elite.front().modal_data;
        }
        saveBestResult();
    }
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total).count();
    
    std::cout << "\n============================================================\n  OPTIMIZATION COMPLETED  \n";
    std::cout << "============================================================\n";
    std::cout << std::fixed << std::setprecision(6) << "  Best error: " << best_error_ << "\n";
    std::cout << std::fixed << std::setprecision(1) << "  Time: " << duration / 1000.0 << " s\n";
    std::cout << "  Modal solver runs: " << evo_params_.population_size * (evo_params_.generations + 1) << "\n";
    std::cout << "============================================================\n";
    
    printFinal(elite.front());
}

// ============================================================
// ВЫВОД РЕЗУЛЬТАТОВ
// ============================================================
void EvolutionaryOptimizer::printHeader() {
    std::cout << "\n";
    std::cout << std::right
              << std::setw(4) << "Gen" << " " << std::setw(12) << "Error"
              << " " << std::setw(10) << "H" << " " << std::setw(10) << "Rm"
              << " " << std::setw(10) << "t_top" << " " << std::setw(10) << "t_bot"
              << " " << std::setw(12) << "Hum" << " " << std::setw(12) << "Prime"
              << " " << std::setw(12) << "Tertia" << " " << std::setw(12) << "Quinta"
              << " " << std::setw(12) << "Nominal" << "\n";
    std::cout << std::string(120, '-') << "\n";
}

void EvolutionaryOptimizer::printGeneration(int gen, const Candidate& best) {
    if (best.error >= 1e8) {
        std::cout << std::right << std::setw(4) << gen << " " << std::setw(12) << "FAILED" << "\n";
        return;
    }
    
    const auto& modes = best.modal_data.modes;
    
    double hum_freq = 0.0, prime_freq = 0.0, tertia_freq = 0.0, quinta_freq = 0.0, nominal_freq = 0.0;
    
    if (modes.size() > 0) hum_freq = modes[0].frequency_hz;
    if (modes.size() > 1) prime_freq = modes[1].frequency_hz;
    if (modes.size() > 2) tertia_freq = modes[2].frequency_hz;
    if (modes.size() > 3) quinta_freq = modes[3].frequency_hz;
    if (modes.size() > 4) nominal_freq = modes[4].frequency_hz;
    
    std::cout << std::right << std::setw(4) << gen
              << std::setw(12) << std::fixed << std::setprecision(4) << best.error
              << std::setw(10) << std::setprecision(1) << best.geometry.H
              << std::setw(10) << best.geometry.R_mouth
              << std::setw(10) << best.geometry.t_top
              << std::setw(10) << best.geometry.t_bottom
              << std::setw(12) << hum_freq
              << std::setw(12) << prime_freq
              << std::setw(12) << tertia_freq
              << std::setw(12) << quinta_freq
              << std::setw(12) << nominal_freq << "\n";
}

void EvolutionaryOptimizer::printFinal(const Candidate& best) {
    if (best.error >= 1e8) {
        std::cout << "\n  No valid solution found.\n";
        return;
    }
    
    std::cout << "\n  Best geometry:\n";
    std::cout << std::fixed << std::setprecision(1) << "    H         = " << best.geometry.H << " mm\n";
    std::cout << std::fixed << std::setprecision(1) << "    R_mouth   = " << best.geometry.R_mouth << " mm\n";
    std::cout << std::fixed << std::setprecision(1) << "    t_top     = " << best.geometry.t_top << " mm\n";
    std::cout << std::fixed << std::setprecision(1) << "    t_bottom  = " << best.geometry.t_bottom << " mm\n";
    std::cout << "    Radii     = ";
    for (double r : best.geometry.r_knots) {
        std::cout << std::fixed << std::setprecision(3) << r << " ";
    }
    std::cout << "\n";
    
    if (best.modal_data.modes.size() >= 5) {
        std::cout << "\n  Frequencies:\n";
        const char* names[] = {"Hum", "Prime", "Tertia", "Quinta", "Nominal"};
        for (int i = 0; i < 5 && i < (int)best.modal_data.modes.size(); ++i) {
            std::cout << std::fixed << std::setprecision(2) << "    " << names[i] << ": " << best.modal_data.modes[i].frequency_hz << " Hz";
            if (i > 0) {
                double ratio = best.modal_data.modes[i].frequency_hz / best.modal_data.modes[1].frequency_hz;
                std::cout << std::fixed << std::setprecision(3) << " (ratio: " << ratio << ")";
            }
            std::cout << "\n";
        }
    }
}