// main.cpp
#include "bell_optimizer.h"
#include "bell_generator.h"
#include "ccx/ccx_solver.hpp"
#include "fem/internal_modal_solver.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>

#include "simple_json.hpp"
#include "models/bell_model_profile.hpp"
#include "bell_material.h"



int command_optimize(int argc, char** argv) {
    std::string workdir;
    std::string evolution_file = "evolution.json";
    std::string targets_file = "targets.json";
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: bell optimize [--solver=ccx|internal] [--workdir <dir>] [--evolution <file>] [--targets <file>]\n";
            return 0;
        }
        if (arg == "--workdir") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --workdir requires a path\n";
                return 2;
            }
            workdir = argv[++i];
            continue;
        }
        if (arg.rfind("--workdir=", 0) == 0) {
            workdir = arg.substr(10);
            if (workdir.empty()) {
                std::cerr << "Error: --workdir requires a non-empty path\n";
                return 2;
            }
            continue;
        }
        if (arg == "--evolution" || arg == "--targets") {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << arg << " requires a path\n";
                return 2;
            }
            if (arg == "--evolution") evolution_file = argv[++i];
            else targets_file = argv[++i];
            continue;
        }
        if (arg.rfind("--evolution=", 0) == 0) {
            evolution_file = arg.substr(12);
            if (evolution_file.empty()) return 2;
            continue;
        }
        if (arg.rfind("--targets=", 0) == 0) {
            targets_file = arg.substr(10);
            if (targets_file.empty()) return 2;
            continue;
        }
        if (arg != "--solver=internal" && arg != "--internal" &&
            arg != "--solver=ccx" && arg != "--ccx") {
            std::cerr << "Error: unknown optimize option: " << arg << "\n";
            return 2;
        }
    }

    if (!workdir.empty()) {
        std::error_code ec;
        std::filesystem::current_path(workdir, ec);
        if (ec) {
            std::cerr << "Error: cannot change working directory to '" << workdir
                      << "': " << ec.message() << "\n";
            return 2;
        }
    }

    // ============================================================
    // 1. МАТЕРИАЛ
    // ============================================================
    const auto& default_material = BellOptimizer::default_material_preset();
    std::string material_id = default_material.id;
    BellMaterial material = default_material.material;

    // ============================================================
    // 2. ЧТЕНИЕ ЦЕЛЕЙ ИЗ JSON
    // ============================================================
    TargetData targets;
    BellOptimizer::Models::BellModelSelection model_selection;
    
    try {
        const auto data = bell::json::parse_file(targets_file);
        model_selection.tuning_model = bell::json::string(data, "tuning_model", "western_five_partial");
        model_selection.family = bell::json::string(data, "tradition", "neutral");
        targets.prime_freq_target = bell::json::number(data, "prime_freq_target", 432.0);
        bell::json::number_array(data, "freq_ratios", targets.freq_ratios.data(), targets.freq_ratios.size());
        bell::json::number_array(data, "amp_ratios", targets.amp_ratios.data(), targets.amp_ratios.size());
        bell::json::number_array(data, "freq_weights", targets.freq_weights.data(), targets.freq_weights.size());
        bell::json::number_array(data, "amp_weights", targets.amp_weights.data(), targets.amp_weights.size());
        std::cout << "✅ Загружены цели из " << targets_file << "\n";
        if (!BellOptimizer::Models::is_supported_tuning_model(model_selection.tuning_model)) {
            std::cerr << "❌ Tuning model '" << model_selection.tuning_model
                      << "' is not supported by the current axisymmetric modal model.\n";
            std::cerr << "   Use western_five_partial, or wait for the non-axisymmetric dual-tone Chinese-bell model.\n";
            return 2;
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Ошибка чтения targets.json: " << e.what() << "\n";
        std::cout << "   Используются значения по умолчанию\n";
    }

    // ============================================================
    // 3. ЧТЕНИЕ ПАРАМЕТРОВ ЭВОЛЮЦИИ ИЗ JSON
    // ============================================================
    EvolutionParams evo;
    
    try {
        const auto data = bell::json::parse_file(evolution_file);
        evo.population_size = bell::json::integer(data, "population_size", 30);
        evo.generations = bell::json::integer(data, "generations", 30);
        evo.num_elite = bell::json::integer(data, "num_elite", 6);
        evo.mutation_probability = bell::json::number(data, "mutation_probability", 0.1);
        evo.mutation_strength = bell::json::number(data, "mutation_strength", 0.05);
        evo.num_knots = bell::json::integer(data, "num_knots", 7);
        evo.seed = static_cast<std::uint64_t>(bell::json::integer(data, "seed", 0));
        evo.deterministic = bell::json::boolean(data, "deterministic", false);
        evo.optimize_ellipticity = bell::json::boolean(data, "optimize_ellipticity", false);
        evo.use_direct_elliptical_fem = bell::json::boolean(data, "use_direct_elliptical_fem", false);
        evo.enable_round_fem_diagnostics = bell::json::boolean(data, "enable_round_fem_diagnostics", false);
        evo.round_fem_diagnostics_csv = bell::json::string(data, "round_fem_diagnostics_csv", "round_fem_diagnostics.csv");
        evo.ellipticity_min = bell::json::number(data, "ellipticity_min", 0.70);
        evo.ellipticity_max = bell::json::number(data, "ellipticity_max", 1.00);
        material_id = bell::json::string(data, "material", material_id);
        const auto* material_preset = BellOptimizer::find_material_preset(material_id);
        if (!material_preset) {
            std::cerr << "Unknown material preset '" << material_id << "'. Available presets:";
            for (const auto& id : BellOptimizer::material_preset_ids()) std::cerr << " " << id;
            std::cerr << "\n";
            return 2;
        }
        material = material_preset->material;
        if (!std::isfinite(evo.ellipticity_min) || !std::isfinite(evo.ellipticity_max) ||
            evo.ellipticity_min < 0.70 || evo.ellipticity_max > 1.0 ||
            evo.ellipticity_min >= evo.ellipticity_max) {
            std::cerr << "Invalid ellipticity search range; require 0.70 <= min < max <= 1.00\n";
            return 2;
        }
        model_selection.family = bell::json::string(data, "tradition", bell::json::string(data, "seed_family", model_selection.family));
        std::cout << "✅ Загружены параметры эволюции из " << evolution_file << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Ошибка чтения evolution.json: " << e.what() << "\n";
        std::cout << "   Используются значения по умолчанию\n";
    }

    // ============================================================
    // 4. ЗАПУСК ОПТИМИЗАЦИИ
    // ============================================================
    std::cout << "\n============================================================\n";
    std::cout << "  BELL OPTIMIZER - Bell Tuning  \n";
    std::cout << "============================================================\n\n";
    
    std::cout << "  Material preset: " << material_id << "\n";
    std::cout << std::fixed << std::setprecision(1) << "  Material: E=" << material.E
              << " MPa, nu=" << std::setprecision(2) << material.nu
              << ", rho=" << std::scientific << std::setprecision(2) << material.rho
              << " t/mm^3\n" << std::defaultfloat;
    std::cout << std::fixed << std::setprecision(1) << "  Target Prime: " << targets.prime_freq_target << " Hz\n";
    std::cout << "  Target ratios (relative to Prime):\n";
    std::cout << std::fixed << std::setprecision(2) << "    Hum     = " << targets.freq_ratios[0] << "  (" << std::setprecision(1) << targets.prime_freq_target * targets.freq_ratios[0] << " Hz)\n";
    std::cout << std::fixed << std::setprecision(2) << "    Tertia  = " << targets.freq_ratios[1] << "  (" << std::setprecision(1) << targets.prime_freq_target * targets.freq_ratios[1] << " Hz)\n";
    std::cout << std::fixed << std::setprecision(2) << "    Quinta  = " << targets.freq_ratios[2] << "  (" << std::setprecision(1) << targets.prime_freq_target * targets.freq_ratios[2] << " Hz)\n";
    std::cout << std::fixed << std::setprecision(2) << "    Nominal = " << targets.freq_ratios[3] << "  (" << std::setprecision(1) << targets.prime_freq_target * targets.freq_ratios[3] << " Hz)\n";
    std::cout << "  Population size: " << evo.population_size << "\n";
    std::cout << "  Generations: " << evo.generations << "\n";
    std::cout << "  Elite count: " << evo.num_elite << "\n";
    std::cout << std::fixed << std::setprecision(0) << "  Mutation probability: " << evo.mutation_probability * 100 << "%\n";
    std::cout << "  Knots: " << evo.num_knots << "\n";
    std::cout << "  Ellipticity search: " << (evo.optimize_ellipticity ? "enabled" : "frozen q=1")
              << " [" << evo.ellipticity_min << ", " << evo.ellipticity_max << "]\n";
    std::cout << "  Direct elliptical vector FEM: "
              << (evo.use_direct_elliptical_fem ? "enabled" : "disabled (legacy morph fallback)") << "\n";
    std::cout << "  RNG: " << (evo.deterministic ? "deterministic" : "random") << ", seed=" << evo.seed << "\n";
    std::cout << "  Model family: " << model_selection.family << "\n";
    std::cout << "  Tuning model: " << model_selection.tuning_model << "\n";
    std::cout << "============================================================\n\n";

    auto start = std::chrono::high_resolution_clock::now();

    // The internal solver is the production default and requires no external
    // executable. CalculiX remains an explicit reference/compatibility path.
    std::string solver_name = "internal";
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--solver=internal" || arg == "--internal") solver_name = "internal";
        else if (arg == "--solver=ccx" || arg == "--ccx") solver_name = "ccx";
    }

    std::unique_ptr<BellOptimizer::ModalSolver> solver;
    if (solver_name == "internal")
        solver = std::make_unique<BellOptimizer::InternalModalSolver>(true, 10, true);
    else
        solver = std::make_unique<BellOptimizer::CalculiXSolver>();
    std::cout << "  Modal solver: " << solver->name() << "\n";
    if (evo.optimize_ellipticity && solver_name != "internal") {
        std::cerr << "Ellipticity optimization requires the internal solver with captured modal displacement fields.\n";
        return 2;
    }

    auto seed_provider = BellOptimizer::Models::make_seed_provider(model_selection.family);
    EvolutionaryOptimizer optimizer(material, targets, evo, std::move(solver), std::move(seed_provider));
    
    // ============================================================
    // 5. ПРОВЕРКА НАЛИЧИЯ СОХРАНЁННОГО РЕЗУЛЬТАТА
    // ============================================================
    std::ifstream check_file("best_geometry.json");
    if (check_file.good()) {
        check_file.close();
        std::cout << "\n📂 Найден сохранённый результат best_geometry.json\n";
        std::cout << "   Загружаем как начальное приближение...\n";
        
        if (optimizer.loadBestFromFile("best_geometry.json")) {
            std::cout << "   ✅ Продолжаем оптимизацию с сохранённого результата\n";
            std::cout << std::fixed << std::setprecision(6) << "   Загружена ошибка: " << optimizer.getBestError() << "\n";
        } else {
            std::cout << "   ⚠️ Не удалось загрузить, начинаем с нуля\n";
        }
    } else {
        std::cout << "\n🆕 Сохранённый результат не найден, начинаем с нуля\n";
    }

    optimizer.run();

    const auto assessment = BellOptimizer::Models::assess_tradition(optimizer.getBestGeometry(), model_selection.family);
    std::cout << "\n  Model geometry assessment [" << model_selection.family << "]:\n";
    std::cout << std::fixed << std::setprecision(1)
              << "    score       = " << assessment.geometry_score << "/100\n"
              << "    proportions = " << assessment.proportion_score * 100.0 << "/100\n"
              << "    thickness   = " << assessment.thickness_score * 100.0 << "/100\n"
              << "    profile     = " << assessment.profile_score * 100.0 << "/100\n"
              << "    verdict     = " << assessment.verdict << "\n";
    if (!assessment.compatible_with_axisymmetric_model)
        std::cout << "    limitation  = " << assessment.limitation << "\n";

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n[TIME] Total: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n\n";

    return 0;
}