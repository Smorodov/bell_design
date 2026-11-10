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
#include "target_profile_loader.hpp"
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
        targets = bell::profiles::load(bell::json::parse_file(targets_file));
        std::cout << "Loaded target profile '" << targets.profile_id << "' (" << targets.display_name << ") from " << targets_file << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Target profile error: " << e.what() << "\n";
        return 2;
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
        const auto* seed_family_value = data.find("seed_family");
        const auto* evolution_tradition_value = data.find("tradition");
        const std::string evolution_tradition = evolution_tradition_value
            ? evolution_tradition_value->string("") : std::string{};
        const std::string configured_seed_family = seed_family_value
            ? seed_family_value->string("")
            : (!evolution_tradition.empty() ? evolution_tradition
               : (!targets.legacy_seed_family.empty() ? targets.legacy_seed_family : model_selection.family));
        if (configured_seed_family.empty()) throw std::runtime_error("seed_family must not be empty");
        if (!targets.legacy_seed_family.empty() &&
            (seed_family_value || evolution_tradition_value) &&
            targets.legacy_seed_family != configured_seed_family)
            throw std::runtime_error("conflicting legacy targets.tradition and evolution geometry-family values");
        if (seed_family_value && !evolution_tradition.empty() &&
            seed_family_value->string("") != evolution_tradition)
            throw std::runtime_error("conflicting evolution.seed_family and legacy evolution.tradition values");
        model_selection.family = configured_seed_family;
        std::cout << "Loaded evolution parameters from " << evolution_file << "\n";
    } catch (const std::exception& e) {
        std::cerr << "❌ Ошибка чтения evolution.json: " << e.what() << "\n";
        return 2;
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
    std::cout << std::fixed << std::setprecision(1)
              << "  Frequency anchor: sorted mode " << targets.frequency_reference_mode_index
              << " -> " << targets.prime_freq_target << " Hz\n";
    std::cout << "  Frequency targets: " << targets.frequency_targets.size()
              << "; amplitude targets: " << targets.amplitude_targets.size() << "\n";
    for (const auto& target : targets.frequency_targets) {
        std::cout << "    mode[" << target.mode_index << "] / reference = " << target.ratio
                  << " (weight " << target.weight << ")\n";
    }
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