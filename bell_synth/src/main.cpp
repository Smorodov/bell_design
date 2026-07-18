// =============================================================================
// main.cpp - СИНТЕЗ КОЛОКОЛА ИЗ .dat ФАЙЛА (С FFT)
// =============================================================================
#include "bell_synthesis.hpp"
#include "../../src/ccx_parser.hpp"

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#ifdef _WIN32
    #include <windows.h>
#endif

using namespace BellOptimizer;
using json = nlohmann::json;

int main(int argc, char* argv[]) {
    #ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
    #endif

    fmt::print("\n{:=^70}\n", "");
    fmt::print("{:=^70}\n", "🔔 СИНТЕЗ КОЛОКОЛА ИЗ .dat (С FFT)");
    fmt::print("{:=^70}\n\n", "");

    std::string filename;
    if (argc > 1) {
        filename = argv[1];
    } else {
        fmt::print("   Введите путь к .dat файлу: ");
        std::getline(std::cin, filename);
        if (filename.empty()) {
            filename = "bell.dat";
        }
    }

    CCXData data = CCXParser::parse_file(filename);
    
    if (data.modes.empty()) {
        fmt::print(stderr, "❌ Ошибка: не найдены моды в файле {}\n", filename);
        return 1;
    }

    // ============================================================
    // ЗАГРУЗКА КОНФИГУРАЦИИ ИЗ JSON
    // ============================================================
    SynthesisConfig config;
    
    try {
        std::ifstream config_file("synthesis.json");
        if (config_file.is_open()) {
            json cfg = json::parse(config_file);
            
            config.sample_rate = cfg.value("sample_rate", 96000);
            config.oversample_factor = cfg.value("oversample_factor", 8);
            config.duration = cfg.value("duration", 0.0);
            config.fade_in = cfg.value("fade_in", 0.01);
            config.fade_out = cfg.value("fade_out", 0.2);
            config.loss_factor = cfg.value("loss_factor", 0.0002);
            config.nonlinear_damping = cfg.value("nonlinear_damping", 0.5);
            config.tail_factor = cfg.value("tail_factor", 3.0);
            config.envelope_threshold = cfg.value("envelope_threshold", 0.001);
            config.boost_amplitude = cfg.value("boost_amplitude", 0.0);
            config.boost_freq_factor = cfg.value("boost_freq_factor", 0.5);
            config.slowdown_factor = cfg.value("slowdown_factor", 0.15);
            config.impact_position = cfg.value("impact_position", 0.70);
            config.impact_duration = cfg.value("impact_duration", 0.0005);
            config.impact_force = cfg.value("impact_force", 1.0);
            config.max_frequency = cfg.value("max_frequency", 20000.0);
            
            config.use_acoustic_model = cfg.value("use_acoustic_model", false);
            config.use_acoustic_coupling = cfg.value("use_acoustic_coupling", false);
            config.use_helmholtz_resonance = cfg.value("use_helmholtz_resonance", false);
            config.helmholtz_amplitude = cfg.value("helmholtz_amplitude", 0.04);
            config.helmholtz_damping = cfg.value("helmholtz_damping", 0.3);
            config.radiation_efficiency = cfg.value("radiation_efficiency", 0.4);
            config.radiation_damping_factor = cfg.value("radiation_damping_factor", 0.005);
            
            config.use_mount_model = cfg.value("use_mount_model", false);
            config.mount_damping = cfg.value("mount_damping", 0.2);
            
            config.reverb_mix = cfg.value("reverb_mix", 0.0);
            config.reverb_delay = cfg.value("reverb_delay", 0.1);
            config.reverb_feedback = cfg.value("reverb_feedback", 0.5);
            
            config.use_rms_normalization = cfg.value("use_rms_normalization", false);
            config.target_rms = cfg.value("target_rms", 0.35);
            config.target_peak = cfg.value("target_peak", 0.95);
            
            fmt::print("✅ Загружена конфигурация из synthesis.json\n");
        } else {
            fmt::print("⚠️ Файл synthesis.json не найден, используются значения по умолчанию\n");
        }
    } catch (const std::exception& e) {
        fmt::print(stderr, "❌ Ошибка чтения synthesis.json: {}\n", e.what());
        fmt::print("   Используются значения по умолчанию\n");
    }

    double impact_force = config.impact_force;
    
    double prime_frequency = 0.0;
    if (data.modes.size() > 1) {
        prime_frequency = data.modes[1].frequency_hz;
    } else {
        prime_frequency = data.modes[0].frequency_hz;
    }
    config.target_fundamental = prime_frequency;
    
    config.alpha_damping = 0.0;
    config.beta_damping = 0.0;

    std::vector<Mode> modes;
    std::vector<double> frequencies;
    std::vector<std::string> mode_names;
    
    double max_amp = 0.0;
    for (const auto& d : data.modes) {
        double amp = std::sqrt(d.effective_mass_norm);
        if (amp > max_amp) max_amp = amp;
    }
    if (max_amp < 1e-10) max_amp = 1.0;
    
    for (size_t i = 0; i < data.modes.size(); ++i) {
        Mode m;
        m.frequency = data.modes[i].frequency_hz;
        frequencies.push_back(m.frequency);
        
        std::string name;
        if (data.modes[i].mode_number == 2) {
            name = "Hum";
        } else if (data.modes[i].mode_number == 3) {
            name = "Prime";
        } else if (data.modes[i].mode_number == 4) {
            name = "Tertia";
        } else if (data.modes[i].mode_number == 5) {
            name = "Quinta";
        } else if (data.modes[i].mode_number == 6) {
            name = "Nominal";
        } else {
            name = "Overtone";
        }
        mode_names.push_back(name);
        
        m.amplitude = std::sqrt(data.modes[i].effective_mass_norm) / max_amp * impact_force;
        
        if (m.frequency > 1.0) {
            m.decay_time = 2.2 / (M_PI * m.frequency * config.loss_factor);
        } else {
            m.decay_time = 5.0;
        }
        
        double freq_ratio_attack = m.frequency / config.target_fundamental;
        m.attack_time = 0.003 + 0.008 / std::sqrt(freq_ratio_attack);
        
        m.decay_time *= impact_force;
        m.attack_time *= impact_force;
        
        if (m.attack_time < 0.001) m.attack_time = 0.001;
        if (m.attack_time > 0.04) m.attack_time = 0.04;
        if (m.decay_time > 30.0) m.decay_time = 30.0;
        if (m.decay_time < 0.1) m.decay_time = 0.1;
        
        m.phase = 0.0;
        modes.push_back(m);
    }

    if (modes.empty()) {
        fmt::print(stderr, "❌ Ошибка: нет мод для синтеза\n");
        return 1;
    }

    // Длительность
    if (config.duration > 0.0) {
        // Используем заданную длительность
        if (config.duration < 0.1) config.duration = 0.1;
        if (config.duration > 30.0) config.duration = 30.0;
    } else {
        // Автоматическая (естественная) длительность
        double max_t60 = 0.0;
        for (const auto& mode : modes) {
            if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
        }
        config.duration = max_t60 * config.tail_factor;
        if (config.duration > 30.0) config.duration = 30.0;
        if (config.duration < 1.0) config.duration = 3.0;
    }

    GeometryParams geo;
    geo.height = 0.45;
    geo.r_max = 0.32;
    geo.r_min = 0.14;
    geo.thickness = 0.022;

    fmt::print("\n📋 ПАРАМЕТРЫ:\n");
    fmt::print("   Основная частота (Prime): {:.1f} Гц\n", config.target_fundamental);
    fmt::print("   Количество мод: {}\n", modes.size());
    fmt::print("   Сила удара: {:.1f}\n", impact_force);
    fmt::print("   Коэффициент потерь: {:.2e}\n", config.loss_factor);
    fmt::print("   Длительность: {:.1f} с\n", config.duration);
    fmt::print("   Fade in: {:.3f} с, Fade out: {:.3f} с\n", config.fade_in, config.fade_out);

    fmt::print("\n🎵 СПЕКТР МОД (ТОЛЬКО ИЗ .dat):\n");
    fmt::print("   {:-<70}\n", "");
    fmt::print("   № | Частота(Гц) | Амплитуда | T60(с) | Название\n");
    fmt::print("   {:-<70}\n", "");
    for (size_t i = 0; i < modes.size() && i < 20; ++i) {
        fmt::print("   {:2d} | {:10.2f} | {:8.6f} | {:6.2f}  | {}\n",
                   i+1, modes[i].frequency, modes[i].amplitude, 
                   modes[i].decay_time, mode_names[i]);
    }
    if (modes.size() > 20) {
        fmt::print("   ... и ещё {} мод\n", modes.size() - 20);
    }
    fmt::print("   {:-<70}\n", "");

    auto start = std::chrono::high_resolution_clock::now();
    
    config.output_filename = fmt::format("bell_{}hz.wav", (int)config.target_fundamental);
    
    BellSynthesizer synth(config);
    synth.generate_modes(geo, modes);
    synth.synthesize();
    
    // Применяем fade из конфигурации
    synth.apply_fade(config.fade_in, config.fade_out);
    
    synth.normalize(config.target_peak);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    fmt::print("\n   Время синтеза: {} мс\n", ms.count());
    synth.save_wav(config.output_filename);
    synth.analyze_spectrum(frequencies);
    synth.print_spectrum(modes);
    
    synth.save_all_envelopes_csv("envelopes_all.csv");
    
    fmt::print("\n   Выполняется FFT анализ...\n");
    
    size_t fft_size = synth.get_audio().size();
    synth.analyze_fft(fft_size);
    
    synth.print_fft_peaks(frequencies);
    synth.save_spectrum_csv(fmt::format("spectrum_{}hz.csv", (int)config.target_fundamental), frequencies);
    
    fmt::print("\n{:=^70}\n", "");
    fmt::print("✅ ГОТОВО: {}\n", config.output_filename);
    fmt::print("   Источник: {}\n", filename);
    fmt::print("   Мод: {}\n", modes.size());
    fmt::print("   Длительность: {:.1f} с\n", synth.get_actual_duration());
    fmt::print("   CSV спектр: spectrum_{}hz.csv\n", (int)config.target_fundamental);
    fmt::print("   CSV огибающих: envelopes_all.csv\n");
    fmt::print("{:=^70}\n\n", "");
    
    return 0;
}