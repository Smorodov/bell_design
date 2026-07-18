// main.cpp
#include "bell_optimizer.h"
#include "bell_generator.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#ifdef _WIN32
#include <windows.h>
#endif

using json = nlohmann::json;

int main() {
    #ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
    #endif

    // ============================================================
    // 1. МАТЕРИАЛ
    // ============================================================
    BellMaterial material;
    material.E = 1.1e4;
    material.nu = 0.34;
    material.rho = 8.6e-9;

    // ============================================================
    // 2. ЧТЕНИЕ ЦЕЛЕЙ ИЗ JSON
    // ============================================================
    TargetData targets;
    
    try {
        std::ifstream targets_file("targets.json");
        if (targets_file.is_open()) {
            json data = json::parse(targets_file);
            
            targets.prime_freq_target = data.value("prime_freq_target", 432.0);
            
            if (data.contains("freq_ratios") && data["freq_ratios"].is_array()) {
                auto arr = data["freq_ratios"];
                for (size_t i = 0; i < 4 && i < arr.size(); ++i) {
                    targets.freq_ratios[i] = arr[i].get<double>();
                }
            }
            
            if (data.contains("amp_ratios") && data["amp_ratios"].is_array()) {
                auto arr = data["amp_ratios"];
                for (size_t i = 0; i < 4 && i < arr.size(); ++i) {
                    targets.amp_ratios[i] = arr[i].get<double>();
                }
            }
            
            if (data.contains("freq_weights") && data["freq_weights"].is_array()) {
                auto arr = data["freq_weights"];
                for (size_t i = 0; i < 4 && i < arr.size(); ++i) {
                    targets.freq_weights[i] = arr[i].get<double>();
                }
            }
            
            if (data.contains("amp_weights") && data["amp_weights"].is_array()) {
                auto arr = data["amp_weights"];
                for (size_t i = 0; i < 4 && i < arr.size(); ++i) {
                    targets.amp_weights[i] = arr[i].get<double>();
                }
            }
            
            targets.prime_weight = data.value("prime_weight", 0.0);
            targets.k_prime = data.value("k_prime", 20.0);
            targets.k_freq = data.value("k_freq", 20.0);
            targets.k_amp = data.value("k_amp", 5.0);
            
            fmt::print("✅ Загружены цели из targets.json\n");
        } else {
            fmt::print("⚠️ Файл targets.json не найден, используются значения по умолчанию\n");
        }
    } catch (const std::exception& e) {
        fmt::print(stderr, "❌ Ошибка чтения targets.json: {}\n", e.what());
        fmt::print("   Используются значения по умолчанию\n");
    }

    // ============================================================
    // 3. ЧТЕНИЕ ПАРАМЕТРОВ ЭВОЛЮЦИИ ИЗ JSON
    // ============================================================
    EvolutionParams evo;
    
    try {
        std::ifstream evo_file("evolution.json");
        if (evo_file.is_open()) {
            json data = json::parse(evo_file);
            
            evo.population_size = data.value("population_size", 30);
            evo.generations = data.value("generations", 30);
            evo.num_elite = data.value("num_elite", 6);
            evo.mutation_probability = data.value("mutation_probability", 0.1);
            evo.mutation_strength = data.value("mutation_strength", 0.05);
            evo.num_knots = data.value("num_knots", 5);
            
            fmt::print("✅ Загружены параметры эволюции из evolution.json\n");
        } else {
            fmt::print("⚠️ Файл evolution.json не найден, используются значения по умолчанию\n");
        }
    } catch (const std::exception& e) {
        fmt::print(stderr, "❌ Ошибка чтения evolution.json: {}\n", e.what());
        fmt::print("   Используются значения по умолчанию\n");
    }

    // ============================================================
    // 4. ЗАПУСК ОПТИМИЗАЦИИ
    // ============================================================
    fmt::print("\n{:=^60}\n", "");
    fmt::print("{:=^60}\n", "  BELL OPTIMIZER - Orthodox Bell Tuning  ");
    fmt::print("{:=^60}\n\n", "");
    
    fmt::print("  Material: E={:.1f} MPa, nu={:.2f}, rho={:.2e} t/mm^3\n", 
               material.E, material.nu, material.rho);
    fmt::print("  Target Prime: {:.1f} Hz\n", targets.prime_freq_target);
    fmt::print("  Target ratios (relative to Prime):\n");
    fmt::print("    Hum     = {:.2f}  ({:.1f} Hz)\n", targets.freq_ratios[0], targets.prime_freq_target * targets.freq_ratios[0]);
    fmt::print("    Tertia  = {:.2f}  ({:.1f} Hz)\n", targets.freq_ratios[1], targets.prime_freq_target * targets.freq_ratios[1]);
    fmt::print("    Quinta  = {:.2f}  ({:.1f} Hz)\n", targets.freq_ratios[2], targets.prime_freq_target * targets.freq_ratios[2]);
    fmt::print("    Nominal = {:.2f}  ({:.1f} Hz)\n", targets.freq_ratios[3], targets.prime_freq_target * targets.freq_ratios[3]);
    fmt::print("  Population size: {}\n", evo.population_size);
    fmt::print("  Generations: {}\n", evo.generations);
    fmt::print("  Elite count: {}\n", evo.num_elite);
    fmt::print("  Mutation probability: {:.0f}%\n", evo.mutation_probability * 100);
    fmt::print("  Knots: {}\n", evo.num_knots);
    fmt::print("{:=^60}\n\n", "");

    auto start = std::chrono::high_resolution_clock::now();

    EvolutionaryOptimizer optimizer(material, targets, evo);
    
    // ============================================================
    // 5. ПРОВЕРКА НАЛИЧИЯ СОХРАНЁННОГО РЕЗУЛЬТАТА
    // ============================================================
    std::ifstream check_file("best_geometry.json");
    if (check_file.good()) {
        check_file.close();
        fmt::print("\n📂 Найден сохранённый результат best_geometry.json\n");
        fmt::print("   Загружаем как начальное приближение...\n");
        
        if (optimizer.loadBestFromFile("best_geometry.json")) {
            fmt::print("   ✅ Продолжаем оптимизацию с сохранённого результата\n");
            fmt::print("   Загружена ошибка: {:.6f}\n", optimizer.getBestError());
        } else {
            fmt::print("   ⚠️ Не удалось загрузить, начинаем с нуля\n");
        }
    } else {
        fmt::print("\n🆕 Сохранённый результат не найден, начинаем с нуля\n");
    }

    optimizer.run();

    auto end = std::chrono::high_resolution_clock::now();
    fmt::print("\n[TIME] Total: {} ms\n\n",
               std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    return 0;
}