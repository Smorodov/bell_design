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
#include "../fem/cylindrical_bell_model.hpp"
#include "../excitation/clapper_modal_state_provider.hpp"
#include "../excitation/clapper_acoustic_objective.hpp"
#include "../acoustics/directional_radiation.hpp"
#include "../excitation/localized_strike.hpp"
#include "../../research/round_bell_audit/round_bell_candidate_analysis.hpp"

namespace {
constexpr double kClapperModalDegeneracyTolerance = 5.0e-5;

std::vector<double> groupedModalFrequencies(const ModalData& data) {
    std::vector<double> frequencies;
    std::vector<std::size_t> counts;
    for (const auto& mode : data.modes) {
        if (!frequencies.empty()) {
            const double previous = frequencies.back();
            const double relative_gap = std::abs(mode.frequency_hz - previous) /
                std::max(mode.frequency_hz, previous);
            if (relative_gap <= kClapperModalDegeneracyTolerance) {
                const std::size_t count = counts.back();
                frequencies.back() += (mode.frequency_hz - frequencies.back()) /
                    static_cast<double>(count + 1);
                counts.back() = count + 1;
                continue;
            }
        }
        frequencies.push_back(mode.frequency_hz);
        counts.push_back(1);
    }
    return frequencies;
}
} // namespace

// ============================================================
// СТРУКТУРА ДЛЯ СОХРАНЕНИЯ РЕЗУЛЬТАТА
// ============================================================
struct BestResult {
    double error = 0.0;
    BellGeometry geometry;
    ModalData modal_data;
    bool group_degenerate_modes = false;
    std::size_t frequency_reference_mode_index = 1;
    std::string objective_model = "frequency_profile";
    
    bool save(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return false;

        file << "{\n";
        file << "    \"error\": " << error << ",\n";
        file << "    \"objective_model\": \"" << objective_model << "\",\n";
        std::vector<double> frequencies;
        if (group_degenerate_modes) frequencies = groupedModalFrequencies(modal_data);
        else for (const auto& mode : modal_data.modes) frequencies.push_back(mode.frequency_hz);
        file << "    \"frequency_reference_mode_index\": " << frequency_reference_mode_index << ",\n";
        file << "    \"frequencies\": {\n";
        if (frequencies.size() > 0) file << "        \"hum\": " << frequencies[0] << ",\n";
        if (frequencies.size() > 1) file << "        \"prime\": " << frequencies[1] << ",\n";
        if (frequencies.size() > 2) file << "        \"tertia\": " << frequencies[2] << ",\n";
        if (frequencies.size() > 3) file << "        \"quinta\": " << frequencies[3] << ",\n";
        if (frequencies.size() > 4) file << "        \"nominal\": " << frequencies[4] << "\n";
        file << "    },\n";
        if (modal_data.clapper_objective_evaluated) {
            file << "    \"clapper\": {\n";
            file << "        \"strike_height_fraction\": " << modal_data.clapper_strike_height_fraction << ",\n";
            file << "        \"mass_kg\": " << modal_data.clapper_mass_kg << ",\n";
            file << "        \"moment_of_inertia_kg_m2\": " << modal_data.clapper_moment_of_inertia_kg_m2 << ",\n";
            file << "        \"com_length_m\": " << modal_data.clapper_com_length_m << ",\n";
            file << "        \"tip_length_m\": " << modal_data.clapper_tip_length_m << ",\n";
            file << "        \"tip_radius_m\": " << modal_data.clapper_tip_radius_m << ",\n";
            file << "        \"hertz_stiffness_n_per_m_1_5\": " << modal_data.clapper_hertz_stiffness_n_per_m_1_5 << ",\n";
            file << "        \"initial_angular_velocity_rad_s\": " << modal_data.clapper_initial_angular_velocity_rad_s << "\n";
            file << "    },\n";
        }
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
    if (targets_.frequency_reference_mode_index < first_data.modes.size()) {
        prime_freq = first_data.modes[targets_.frequency_reference_mode_index].frequency_hz;
    }
    
    if (prime_freq < 1.0) {
        std::cout << std::fixed << std::setprecision(2) << "   ❌ Invalid prime frequency: " << prime_freq << "\n";
        return false;
    }
    
    std::cout << std::fixed << std::setprecision(2) << "   Current reference mode[" << targets_.frequency_reference_mode_index << "]: " << prime_freq << " Hz\n";
    
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
    if (evo_params_.objective_model == "clapper_acoustic") {
        best_error_ = evaluateCandidate(scaled_geometry, scaled_data);
    } else {
        best_error_ = computeError(scaled_data, scaled_geometry);
    }
    best_geometry_ = scaled_geometry;
    best_modal_data_ = scaled_data;
    modal_reference_ = scaled_data;
    has_modal_reference_ = !modal_reference_.modes.empty();
    
    std::cout << "✅ Rescaled successfully!\n";
    std::cout << std::fixed << std::setprecision(6) << "   New error: " << best_error_ << "\n";
    if (targets_.frequency_reference_mode_index < best_modal_data_.modes.size())
        std::cout << std::fixed << std::setprecision(2) << "   New reference mode[" << targets_.frequency_reference_mode_index << "]: " << best_modal_data_.modes[targets_.frequency_reference_mode_index].frequency_hz << " Hz\n";
    
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
    result.group_degenerate_modes = evo_params_.objective_model == "clapper_acoustic";
    result.frequency_reference_mode_index = targets_.frequency_reference_mode_index;
    result.objective_model = evo_params_.objective_model;
    
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
    const auto& modes = modal_data.modes;
    const std::size_t reference = targets_.frequency_reference_mode_index;
    if (reference >= modes.size() || modes[reference].frequency_hz < BellConstants::MIN_FREQ_HZ)
        return BellConstants::ERROR_INVALID;
    const double reference_freq = modes[reference].frequency_hz;

    const auto& frequency_targets = targets_.frequency_targets;
    const auto& amplitude_targets = targets_.amplitude_targets;
    if (frequency_targets.empty() && amplitude_targets.empty()) return BellConstants::ERROR_INVALID;
    double freq_error = 0.0;
    for (const auto& target : frequency_targets) {
        if (target.weight == 0.0) continue;
        if (target.mode_index >= modes.size()) return BellConstants::ERROR_INVALID;
        const double actual = modes[target.mode_index].frequency_hz / reference_freq;
        if (!(actual > 0.0) || !std::isfinite(actual)) return BellConstants::ERROR_INVALID;
        const double e = std::log(actual / target.ratio);
        freq_error += target.weight * e * e;
    }

    double amp_error = 0.0;
    for (const auto& target : amplitude_targets) {
        if (target.weight == 0.0) continue;
        if (target.mode_index >= modes.size() || target.reference_mode_index >= modes.size())
            return BellConstants::ERROR_INVALID;
        const double numerator = modes[target.mode_index].amplitude;
        const double denominator = modes[target.reference_mode_index].amplitude;
        if (!(numerator > 1e-12) || !(denominator > 1e-12) ||
            !std::isfinite(numerator) || !std::isfinite(denominator))
            return BellConstants::ERROR_INVALID;
        const double actual = numerator / denominator;
        if (!(actual > 0.0) || !std::isfinite(actual)) return BellConstants::ERROR_INVALID;
        const double e = std::log(actual / target.ratio);
        amp_error += target.weight * e * e;
    }

    // Keep the established smooth proportion penalty independent of profile metadata.
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
    if (evo_params_.objective_model == "clapper_acoustic") geometry.ellipticity = 1.0;
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
    
    if (evo_params_.objective_model == "clapper_acoustic") {
        try {
            BellMeshConfig fem_mesh;
            // This is the inner optimization loop, so use a deliberately coarse
            // but complete circular modal basis. Final verification/rendering uses
            // the higher-resolution settings outside the evolutionary loop.
            fem_mesh.nNodes = 11;
            FEM::CylindricalModalEvaluationConfig fem_config;
            fem_config.harmonics = {0, 1, 2};
            fem_config.angular_samples = 48;
            fem_config.requested_modes = 14;
            fem_config.max_iterations = 700;
            fem_config.tolerance = 1.0e-6;
            fem_config.minimum_frequency_hz = BellConstants::MIN_FREQ_HZ;

            auto solve_round_modes = [&](const BellGeometry& candidate, ModalData& modes,
                                         FEM::CylindricalBellModel& model) -> bool {
                model = FEM::make_cylindrical_bell_model(candidate, material_, fem_mesh);
                ModalInfluence neutral;
                neutral.radial.assign(model.mesh_si.nodes.size(), 0.0);
                neutral.axial.assign(model.mesh_si.nodes.size(), 0.0);
                return FEM::evaluate_cylindrical_bell_modes(model, neutral, fem_config, modes);
            };

            // First obtain the frequency anchor on the canonical shape, grouping
            // numerically degenerate circular-bell pairs exactly as the objective does.
            ModalData anchor_modes;
            FEM::CylindricalBellModel anchor_model;
            if (!solve_round_modes(geometry, anchor_modes, anchor_model)) { std::cerr << "[clapper-objective] anchor FEM solve failed\n"; return BellConstants::ERROR_INVALID; }
            const std::vector<double> grouped_frequencies = groupedModalFrequencies(anchor_modes);
            if (targets_.frequency_reference_mode_index >= grouped_frequencies.size()) {
                std::cerr << "[clapper-objective] frequency reference exceeds grouped mode count\n";
                return BellConstants::ERROR_INVALID;
            }
            const double anchor_frequency = grouped_frequencies[targets_.frequency_reference_mode_index];
            const double physical_scale = anchor_frequency / targets_.prime_freq_target;
            if (!std::isfinite(physical_scale) || physical_scale < 0.05 || physical_scale > 20.0)
                return BellConstants::ERROR_INVALID;
            geometry.scale(physical_scale);

            const double top_ratio = geometry.t_top / geometry.H;
            const double bottom_ratio = geometry.t_bottom / geometry.H;
            const double diameter_ratio = 2.0 * geometry.R_mouth / geometry.H;
            if (top_ratio < 0.005 || top_ratio > 0.050 || bottom_ratio < 0.010 || bottom_ratio > 0.100 ||
                diameter_ratio < 0.7 || diameter_ratio > 1.3) return BellConstants::ERROR_INVALID;

            // The contact model must use FEM modes and the inner profile from this
            // final, physically scaled geometry. No pre-scale boundary data is reused.
            ModalData final_modes;
            FEM::CylindricalBellModel final_model;
            if (!solve_round_modes(geometry, final_modes, final_model)) { std::cerr << "[clapper-objective] final-geometry FEM solve failed\n"; return BellConstants::ERROR_INVALID; }
            if (final_model.mesh_si.inner_profile.size() < 2 || final_model.mesh_si.outer_profile.size() < 2)
                return BellConstants::ERROR_INVALID;

            std::vector<Point2D> outer_profile_mm = final_model.mesh_si.outer_profile;
            for (auto& point : outer_profile_mm) { point.r *= 1000.0; point.z *= 1000.0; }
            for (auto& mode : final_modes.modes) {
                const auto field = Acoustics::buildDirectionalRadiationField(
                    mode, outer_profile_mm, evo_params_.clapper.listener_polar_rad, 1.0);
                mode.acoustic_pressure_coefficient = std::abs(
                    field.evaluate(evo_params_.clapper.listener_azimuth_rad));
            }

            Excitation::ClapperModalEvaluationConfig clapper_eval;
            clapper_eval.top_z_m = final_model.mesh_si.inner_profile.front().z;
            clapper_eval.mouth_z_m = final_model.mesh_si.inner_profile.back().z;
            // The present reference suspension is demonstrative rather than
            // calibrated. Accept a broad physically plausible strike band while
            // recording the actual fraction in the result diagnostics.
            clapper_eval.admissible_strike_fraction_min = 0.35;
            clapper_eval.admissible_strike_fraction_max = 0.70;
            if (!Excitation::makeScaledClapperDynamicsConfig(
                    geometry.H * 1.0e-3, clapper_eval.top_z_m, evo_params_.clapper,
                    material_, clapper_eval.clapper)) { std::cerr << "[clapper-objective] clapper scaling failed\n"; return BellConstants::ERROR_INVALID; }
            Excitation::ClapperModalEvaluationResult excitation;
            if (!Excitation::evaluateClapperModalStateFromCandidate(
                    clapper_eval, final_model.mesh_si.inner_profile, final_modes.modes, excitation)) {
                std::cerr << "[clapper-objective] clapper contact/evaluation failed (contact="
                          << excitation.contact.contact_occurred << ", strike_fraction="
                          << excitation.strike_height_fraction << ")\n";
                return BellConstants::ERROR_INVALID;
            }
            final_modes.clapper_objective_evaluated = true;
            final_modes.clapper_strike_height_fraction = excitation.strike_height_fraction;
            final_modes.clapper_mass_kg = clapper_eval.clapper.mass_kg;
            final_modes.clapper_moment_of_inertia_kg_m2 = clapper_eval.clapper.moment_of_inertia_kg_m2;
            final_modes.clapper_com_length_m = clapper_eval.clapper.centre_of_mass_length_m;
            final_modes.clapper_tip_length_m = clapper_eval.clapper.tip_length_m;
            final_modes.clapper_tip_radius_m = clapper_eval.clapper.tip_radius_m;
            final_modes.clapper_hertz_stiffness_n_per_m_1_5 = clapper_eval.clapper.hertz_stiffness_n_per_m_1_5;
            final_modes.clapper_initial_angular_velocity_rad_s = clapper_eval.clapper.initial_angular_velocity_rad_s;

            Excitation::ClapperModalState state;
            state.modal_data = final_modes;
            state.inner_profile_si = final_model.mesh_si.inner_profile;
            state.excitation = std::move(excitation);
            Excitation::ClapperAcousticObjectiveResult acoustic_score;
            if (!Excitation::evaluateClapperAcousticObjective(state, targets_, acoustic_score)) {
                std::cerr << "[clapper-objective] acoustic target evaluation failed\n";
                return BellConstants::ERROR_INVALID;
            }
            out_data = std::move(final_modes);
            return acoustic_score.total_error;
        } catch (const std::exception& error) {
            std::cerr << "[clapper-objective] exception: " << error.what() << "\n";
            return BellConstants::ERROR_INVALID;
        }
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
    // Optional circular-bell integration side channel. This direct FEM solve
    // is deliberately diagnostic-only: it does not replace raw_data and never
    // contributes to the production objective until modal mapping is validated.
    if (evo_params_.enable_round_fem_diagnostics &&
        std::abs(geometry.ellipticity - 1.0) <= 1.0e-12) {
        try {
            const auto fem_model = FEM::make_cylindrical_bell_model(geometry, material_, config);
            FEM::CylindricalModalEvaluationConfig fem_config;
            fem_config.harmonics = {0, 1, 2, 3};
            fem_config.angular_samples = 64;
            fem_config.requested_modes = 18;
            fem_config.max_iterations = 1600;
            fem_config.tolerance = 1.0e-6;
            fem_config.minimum_frequency_hz = BellConstants::MIN_FREQ_HZ;
            Excitation::MeridionalInfluence fem_influence;
            ModalData diagnostic_modes;
            if (!FEM::evaluate_cylindrical_bell_modes(
                    fem_model, fem_influence.make(fem_model.mesh_si), fem_config, diagnostic_modes)) {
                std::cerr << "[round-fem-diagnostics] direct FEM solve failed; production candidate unchanged\n";
            } else {
                double z_min = fem_model.mesh_si.outer_profile.front().z;
                double z_max = z_min;
                for (const auto& point : fem_model.mesh_si.outer_profile) {
                    z_min = std::min(z_min, point.z);
                    z_max = std::max(z_max, point.z);
                }
                Excitation::CircularLocalizedStrikeConfig strike;
                strike.center_z = z_min + 0.70 * (z_max - z_min);
                strike.meridional_sigma = 0.012;
                strike.circumferential_sigma = 0.008;
                Research::ModalFamilyAdapterConfig adapter_config;
                Research::RoundBellCandidateAnalysis analysis;
                if (!Research::analyzeRoundBellCandidate(
                        fem_model.mesh_si, diagnostic_modes, strike, adapter_config, analysis)) {
                    std::cerr << "[round-fem-diagnostics] modal/strike analysis rejected; production candidate unchanged\n";
                } else {
                    std::ofstream diagnostic_csv(evo_params_.round_fem_diagnostics_csv, std::ios::app);
                    if (diagnostic_csv) {
                        if (diagnostic_csv.tellp() == std::streampos(0))
                            diagnostic_csv << "H_mm,R_mouth_mm,ellipticity,mode_count,family_count,max_abs_strike_projection\n";
                        diagnostic_csv << std::setprecision(12) << geometry.H << ',' << geometry.R_mouth << ','
                                       << geometry.ellipticity << ',' << diagnostic_modes.modes.size() << ','
                                       << analysis.modal_families.size() << ','
                                       << analysis.maximum_absolute_strike_projection << '\n';
                    } else {
                        std::cerr << "[round-fem-diagnostics] cannot open CSV; production candidate unchanged\n";
                    }
                }
            }
        } catch (const std::exception& error) {
            std::cerr << "[round-fem-diagnostics] " << error.what()
                      << "; production candidate unchanged\n";
        }
    }

    // Track modal identity only after the selected geometry path has produced
    // its final modal set. Earlier tracking was silently undone by the direct
    // elliptical FEM replacement and by the legacy morph's frequency sort.
    // The reference is established from the first valid transformed solve and
    // remains fixed during this optimization run.
    if (!has_modal_reference_) {
        modal_reference_ = raw_data;
        has_modal_reference_ = !modal_reference_.modes.empty();
    } else {
        const auto tracking = trackModalFamilies(modal_reference_, raw_data);
        if (tracking.used_shape_tracking && tracking.ambiguous) {
            // Keep a physically valid candidate, but make the assignment
            // uncertainty visible instead of silently treating rank as identity.
            std::cout << "   ⚠️ Modal-family tracking: near-degenerate/ambiguous pair, "
                      << "min shape similarity=" << std::fixed << std::setprecision(3)
                      << tracking.minimum_similarity << "\n" << std::defaultfloat;
        }
    }

    double prime_freq = 0.0;
    if (targets_.frequency_reference_mode_index < raw_data.modes.size()) {
        prime_freq = raw_data.modes[targets_.frequency_reference_mode_index].frequency_hz;
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
    // ЭТАП 8: SELECT THE FINAL-GEOMETRY OBJECTIVE
    // ============================================================
    if (evo_params_.objective_model != "frequency_profile") return BellConstants::ERROR_INVALID;
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

// Accurate second stage for circular bells. It deliberately uses only the
// axisymmetric (m=0) cylindrical FEM family, which is the family audited against
// the fast meridional solver in Stage 13. Non-axisymmetric families are not
// mixed into the frequency-index objective.
double EvolutionaryOptimizer::evaluateRoundCandidateAccurate(
    BellGeometry& geometry, ModalData& out_data) {
    if (evo_params_.objective_model != "frequency_profile" ||
        std::abs(geometry.ellipticity - 1.0) > 1.0e-12)
        return BellConstants::ERROR_INVALID;

    try {
        geometry = canonicalizeShape(geometry);
        BellMeshConfig mesh;
        mesh.nNodes = 15;
        const auto model = FEM::make_cylindrical_bell_model(geometry, material_, mesh);
        FEM::CylindricalModalEvaluationConfig config;
        config.harmonics = {0};
        config.angular_samples = 48;
        config.requested_modes = 18;
        config.max_iterations = 1200;
        config.tolerance = 1.0e-6;
        config.minimum_frequency_hz = BellConstants::MIN_FREQ_HZ;
        Excitation::MeridionalInfluence fem_influence;
        const auto influence = fem_influence.make(model.mesh_si);
        ModalData modes;
        const std::size_t reference = targets_.frequency_reference_mode_index;
        double scale = 0.0;
        bool covers_synthesis_band = false;
        while (true) {
            if (!FEM::evaluate_cylindrical_bell_modes(model, influence, config, modes) || modes.modes.empty())
                return BellConstants::ERROR_INVALID;
            modes.sort_by_frequency();
            if (reference >= modes.modes.size()) return BellConstants::ERROR_INVALID;
            const double prime = modes.modes[reference].frequency_hz;
            if (!std::isfinite(prime) || prime < BellConstants::MIN_FREQ_HZ)
                return BellConstants::ERROR_INVALID;
            scale = prime / targets_.prime_freq_target;
            if (!std::isfinite(scale) || scale < 0.05 || scale > 20.0)
                return BellConstants::ERROR_INVALID;
            const double physical_top = modes.modes.back().frequency_hz / scale;
            covers_synthesis_band = physical_top >= evo_params_.round_synthesis_max_frequency_hz;
            if (covers_synthesis_band || config.requested_modes >= evo_params_.round_synthesis_max_modes ||
                modes.modes.size() < static_cast<std::size_t>(config.requested_modes))
                break;
            config.requested_modes = std::min(evo_params_.round_synthesis_max_modes,
                                               config.requested_modes + 8);
        }
        if (!covers_synthesis_band) {
            std::cout << "[round-modal-band] requested mode cap reached; physical spectrum covered to "
                      << modes.modes.back().frequency_hz / scale << " Hz of "
                      << evo_params_.round_synthesis_max_frequency_hz << " Hz\n";
        }
        geometry.scale(scale);
        for (auto& mode : modes.modes) mode.frequency_hz /= scale;
        const double top = geometry.t_top / geometry.H;
        const double bottom = geometry.t_bottom / geometry.H;
        const double diameter = 2.0 * geometry.R_mouth / geometry.H;
        if (top < 0.005 || top > 0.050 || bottom < 0.010 || bottom > 0.100 ||
            diameter < 0.7 || diameter > 1.3) return BellConstants::ERROR_INVALID;
        const double score = computeError(modes, geometry);
        if (!std::isfinite(score) || score >= 1.0e8) return BellConstants::ERROR_INVALID;
        out_data = std::move(modes);
        return score;
    } catch (const std::exception& error) {
        std::cerr << "[round-two-stage] accurate FEM refinement failed: "
                  << error.what() << "\n";
        return BellConstants::ERROR_INVALID;
    }
}

// ============================================================
// ОЦЕНКА ПОПУЛЯЦИИ: FAST PASS, THEN SELECTIVE FEM REFINEMENT
// ============================================================
void EvolutionaryOptimizer::evaluatePopulation(std::vector<Candidate>& population) {
    const bool two_stage = evo_params_.enable_round_two_stage_optimization &&
        evo_params_.objective_model == "frequency_profile";
    std::vector<std::size_t> evaluated;
    std::vector<BellGeometry> original_geometry(population.size());

    for (std::size_t i = 0; i < population.size(); ++i) {
        auto& cand = population[i];
        // A loaded/best candidate may carry a valid score from another target or
        // an earlier fidelity. Re-evaluate it in this run so it cannot bypass FEM.
        if (!two_stage && cand.valid) continue;
        original_geometry[i] = cand.geometry;
        cand.accurate_evaluated = false;
        ModalData modal_data;
        cand.error = evaluateCandidate(cand.geometry, modal_data);
        cand.modal_data = std::move(modal_data);
        cand.valid = std::isfinite(cand.error) && cand.error < 1.0e8;
        evaluated.push_back(i);
    }

    if (!two_stage || evo_params_.round_two_stage_shortlist_size <= 0 || evaluated.empty())
        return;

    std::vector<std::size_t> ranked;
    ranked.reserve(evaluated.size());
    for (const auto i : evaluated) {
        if (population[i].valid &&
            std::abs(original_geometry[i].ellipticity - 1.0) <= 1.0e-12)
            ranked.push_back(i);
    }
    std::stable_sort(ranked.begin(), ranked.end(), [&](std::size_t a, std::size_t b) {
        if (population[a].error != population[b].error)
            return population[a].error < population[b].error;
        return a < b;
    });

    // At least num_elite candidates must be offered to the accurate stage;
    // otherwise the next generation could have too few accurate parents.
    const std::size_t requested = static_cast<std::size_t>(std::max(
        evo_params_.round_two_stage_shortlist_size, evo_params_.num_elite));
    const std::size_t count = std::min(requested, ranked.size());
    std::size_t refined = 0;
    for (std::size_t n = 0; n < count; ++n) {
        const std::size_t i = ranked[n];
        BellGeometry refined_geometry = original_geometry[i];
        ModalData refined_modes;
        const double refined_error = evaluateRoundCandidateAccurate(refined_geometry, refined_modes);
        if (std::isfinite(refined_error) && refined_error < 1.0e8) {
            population[i].geometry = std::move(refined_geometry);
            population[i].modal_data = std::move(refined_modes);
            population[i].error = refined_error;
            population[i].accurate_evaluated = true;
            ++refined;
        }
    }
    std::cout << "[round-two-stage] fast candidates=" << evaluated.size()
              << ", FEM shortlist=" << count << ", refined=" << refined
              << "; elite selection uses one fidelity only\n";
}

// ============================================================
// ОТБОР
// ============================================================
std::vector<EvolutionaryOptimizer::Candidate> EvolutionaryOptimizer::selectElite(const std::vector<Candidate>& population) {
    std::vector<Candidate> sorted;
    const bool prefer_accurate = evo_params_.enable_round_two_stage_optimization &&
        evo_params_.objective_model == "frequency_profile" &&
        std::any_of(population.begin(), population.end(),
            [](const Candidate& candidate) { return candidate.accurate_evaluated; });
    for (const auto& candidate : population) {
        if (!candidate.valid || !std::isfinite(candidate.error) || candidate.error >= 1.0e8)
            continue;
        if (!prefer_accurate || candidate.accurate_evaluated)
            sorted.push_back(candidate);
    }
    std::stable_sort(sorted.begin(), sorted.end(),
        [](const Candidate& a, const Candidate& b) { return a.error < b.error; });

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
    std::cout << "\n============================================================\n  EVOLUTIONARY OPTIMIZER (Generic target-profile objective)  \n";
    std::cout << "============================================================\n";
    std::cout << std::fixed << std::setprecision(1) << "  Target Prime: " << targets_.prime_freq_target << " Hz\n";
    std::cout << "  Population: " << evo_params_.population_size << "\n";
    std::cout << "  Generations: " << evo_params_.generations << "\n";
    std::cout << std::fixed << std::setprecision(0) << "  Mutation: " << evo_params_.mutation_probability * 100 << "%\n";
    if (evo_params_.objective_model == "clapper_acoustic")
        std::cout << "  Mode: final-geometry cylindrical FEM + scaled clapper contact + acoustic profile score\n";
    else
        std::cout << "  Mode: one axisymmetric modal solve per candidate + analytical frequency scaling\n";
    std::cout << "  Ellipticity search: " << (evo_params_.optimize_ellipticity ? "enabled" : "frozen at q=1") << "\n";
    std::cout << "============================================================\n\n";
    
    auto start_total = std::chrono::high_resolution_clock::now();
    
    std::vector<Candidate> population = initializePopulation();
    evaluatePopulation(population);
    
    auto elite = selectElite(population);
    if (elite.empty()) {
        std::cout << "No valid candidates survived evaluation; optimization stopped safely.\n";
        best_geometry_ = getDefaultGeometry();
        best_error_ = BellConstants::ERROR_INVALID;
        return;
    }
    
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
    std::cout << "  Candidate objective: " << evo_params_.objective_model << "\n";
    std::cout << "============================================================\n";
    
    printFinal(elite.front());
}

// ============================================================
// ВЫВОД РЕЗУЛЬТАТОВ
// ============================================================
void EvolutionaryOptimizer::printHeader() {
    std::cout << "\n" << std::right
              << std::setw(4) << "Gen" << " " << std::setw(12) << "Error"
              << " " << std::setw(10) << "H" << " " << std::setw(10) << "Rm"
              << " " << std::setw(10) << "t_top" << " " << std::setw(10) << "t_bot"
              << " " << std::setw(12) << "Mode0"
              << " " << std::setw(12) << (evo_params_.objective_model == "clapper_acoustic" ? "RefGroup" : "Reference") << "\n";
    std::cout << std::string(100, '-') << "\n";
}

void EvolutionaryOptimizer::printGeneration(int gen, const Candidate& best) {
    if (best.error >= BellConstants::ERROR_THRESHOLD_VALID) {
        std::cout << std::right << std::setw(4) << gen << " " << std::setw(12) << "FAILED" << "\n";
        return;
    }
    const auto& modes = best.modal_data.modes;
    double mode0 = modes.empty() ? 0.0 : modes[0].frequency_hz;
    double reference = targets_.frequency_reference_mode_index < modes.size()
        ? modes[targets_.frequency_reference_mode_index].frequency_hz : 0.0;
    if (evo_params_.objective_model == "clapper_acoustic") {
        const auto groups = groupedModalFrequencies(best.modal_data);
        mode0 = groups.empty() ? 0.0 : groups.front();
        reference = targets_.frequency_reference_mode_index < groups.size()
            ? groups[targets_.frequency_reference_mode_index] : 0.0;
    }
    std::cout << std::right << std::setw(4) << gen
              << std::setw(12) << std::fixed << std::setprecision(4) << best.error
              << std::setw(10) << std::setprecision(1) << best.geometry.H
              << std::setw(10) << best.geometry.R_mouth
              << std::setw(10) << best.geometry.t_top
              << std::setw(10) << best.geometry.t_bottom
              << std::setw(12) << mode0
              << std::setw(12) << reference << "\n";
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
    for (double r : best.geometry.r_knots) std::cout << std::fixed << std::setprecision(3) << r << " ";
    std::cout << "\n";

    if (!best.modal_data.modes.empty()) {
        std::vector<double> displayed_frequencies;
        if (best.modal_data.clapper_objective_evaluated) {
            displayed_frequencies = groupedModalFrequencies(best.modal_data);
            std::cout << "\n  Degeneracy-grouped modal frequencies (objective selector slots):\n";
        } else {
            for (const auto& mode : best.modal_data.modes)
                displayed_frequencies.push_back(mode.frequency_hz);
            std::cout << "\n  Sorted modal frequencies:\n";
        }
        const std::size_t count = std::min<std::size_t>(5, displayed_frequencies.size());
        for (std::size_t i = 0; i < count; ++i) {
            std::cout << std::fixed << std::setprecision(2) << "    mode[" << i << "]: "
                      << displayed_frequencies[i] << " Hz";
            if (i == targets_.frequency_reference_mode_index)
                std::cout << " (frequency reference)";
            std::cout << "\n";
        }
    }

    if (best.modal_data.clapper_objective_evaluated) {
        const auto& m = best.modal_data;
        std::cout << "\n  Geometrically scaled clapper (reference design, not calibrated):\n"
                  << std::fixed << std::setprecision(4)
                  << "    strike fraction = " << m.clapper_strike_height_fraction << "\n"
                  << "    mass            = " << m.clapper_mass_kg << " kg\n"
                  << "    COM length      = " << m.clapper_com_length_m << " m\n"
                  << "    tip length      = " << m.clapper_tip_length_m << " m\n"
                  << "    tip radius      = " << m.clapper_tip_radius_m << " m\n"
                  << "    inertia         = " << m.clapper_moment_of_inertia_kg_m2 << " kg m^2\n"
                  << "    Hertz stiffness = " << m.clapper_hertz_stiffness_n_per_m_1_5 << " N/m^1.5\n"
                  << "    initial omega   = " << m.clapper_initial_angular_velocity_rad_s << " rad/s\n";
    }
}
