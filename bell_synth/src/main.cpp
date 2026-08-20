// =============================================================================
// main.cpp - СИНТЕЗ КОЛОКОЛА (POLYBLEP + REVERB + STEREO)
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
    fmt::print("{:=^70}\n", "🔔 СИНТЕЗ КОЛОКОЛА (POLYBLEP + REVERB + STEREO)");
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
            
            // Основные параметры
            config.sample_rate = cfg.value("sample_rate", 48000);
            config.oversample_factor = cfg.value("oversample_factor", 4);
            config.duration = cfg.value("duration", 0.0);
            config.fade_in = cfg.value("fade_in", 0.003);
            config.fade_out = cfg.value("fade_out", 0.3);
            config.loss_factor = cfg.value("loss_factor", 0.0005);
            config.nonlinear_damping = cfg.value("nonlinear_damping", 0.6);
            config.tail_factor = cfg.value("tail_factor", 1.2);
            config.impact_force = cfg.value("impact_force", 1.0);
            config.target_peak = cfg.value("target_peak", 0.95);
            
            // Параметры обрезки тишины
            config.envelope_threshold = cfg.value("envelope_threshold", 0.00001);
            config.use_tail_trimming = cfg.value("use_tail_trimming", true);
            
            // Реверберация
            config.reverb_mix = cfg.value("reverb_mix", 0.05);
            config.reverb_delay = cfg.value("reverb_delay", 0.02);
            config.reverb_feedback = cfg.value("reverb_feedback", 0.3);
            
            // Стерео-модуляция
            config.use_stereo_modulation = cfg.value("use_stereo_modulation", false);
            config.stereo_depth = cfg.value("stereo_depth", 0.002);
            config.stereo_rate = cfg.value("stereo_rate", 0.15);
            config.stereo_mix = cfg.value("stereo_mix", 0.25);
            
            // Waveform type
            std::string waveform_str = cfg.value("waveform_type", "mixed");
            if (waveform_str == "sine") {
                config.waveform_type = SynthesisConfig::WaveformType::SINE;
            } else if (waveform_str == "triangle") {
                config.waveform_type = SynthesisConfig::WaveformType::TRIANGLE;
            } else if (waveform_str == "square") {
                config.waveform_type = SynthesisConfig::WaveformType::SQUARE;
            } else if (waveform_str == "sawtooth") {
                config.waveform_type = SynthesisConfig::WaveformType::SAWTOOTH;
            } else {
                config.waveform_type = SynthesisConfig::WaveformType::MIXED;
            }
            
            fmt::print("✅ Загружена конфигурация из synthesis.json\n");
            fmt::print("   Waveform: {}\n", waveform_str);
            fmt::print("   Reverb: {:.0f}%\n", config.reverb_mix * 100.0);
            if (config.use_stereo_modulation) {
                fmt::print("   Stereo modulation: ON (depth={:.1f}ms, rate={:.2f}Hz, mix={:.0f}%)\n",
                           config.stereo_depth * 1000.0, config.stereo_rate, config.stereo_mix * 100.0);
            } else {
                fmt::print("   Stereo modulation: OFF\n");
            }
        } else {
            fmt::print("⚠️ Файл synthesis.json не найден, используются значения по умолчанию\n");
        }
    } catch (const std::exception& e) {
        fmt::print(stderr, "❌ Ошибка чтения synthesis.json: {}\n", e.what());
        fmt::print("   Используются значения по умолчанию\n");
    }

    // Определяем основную частоту
    double prime_frequency = 0.0;
    if (data.modes.size() > 1) {
        prime_frequency = data.modes[1].frequency_hz;
    } else {
        prime_frequency = data.modes[0].frequency_hz;
    }
    config.target_fundamental = prime_frequency;

    // ============================================================
    // ПОДГОТОВКА МОД
    // ============================================================
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
        
        m.amplitude = std::sqrt(data.modes[i].effective_mass_norm) / max_amp * config.impact_force;
        
        if (m.frequency > 1.0) {
            m.decay_time = 2.2 / (M_PI * m.frequency * config.loss_factor);
        } else {
            m.decay_time = 5.0;
        }
        
        double freq_ratio_attack = m.frequency / config.target_fundamental;
        m.attack_time = 0.002 + 0.005 / std::sqrt(freq_ratio_attack);
        
        m.decay_time *= config.impact_force;
        m.attack_time *= config.impact_force;
        
        m.attack_time = std::clamp(m.attack_time, 0.001, 0.03);
        m.decay_time = std::clamp(m.decay_time, 0.1, 30.0);
        
        m.phase = 0.0;
        modes.push_back(m);
    }

    if (modes.empty()) {
        fmt::print(stderr, "❌ Ошибка: нет мод для синтеза\n");
        return 1;
    }

    // Длительность
    if (config.duration > 0.0) {
        config.duration = std::clamp(config.duration, 0.1, 30.0);
    } else {
        double max_t60 = 0.0;
        for (const auto& mode : modes) {
            if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
        }
        config.duration = max_t60 * config.tail_factor;
        config.duration = std::clamp(config.duration, 1.0, 30.0);
    }

    GeometryParams geo;
    geo.height = 0.45;
    geo.r_max = 0.32;
    geo.r_min = 0.14;
    geo.thickness = 0.022;

    fmt::print("\n📋 ПАРАМЕТРЫ:\n");
    fmt::print("   Основная частота (Prime): {:.1f} Гц\n", config.target_fundamental);
    fmt::print("   Количество мод: {}\n", modes.size());
    fmt::print("   Сила удара: {:.1f}\n", config.impact_force);
    fmt::print("   Коэффициент потерь: {:.2e}\n", config.loss_factor);
    fmt::print("   Длительность: {:.1f} с\n", config.duration);
    fmt::print("   Реверберация: {:.0f}%\n", config.reverb_mix * 100.0);
    if (config.use_stereo_modulation) {
        fmt::print("   Стерео-модуляция: ВКЛ (глубина={:.1f}мс, частота={:.2f}Гц, уровень={:.0f}%)\n",
                   config.stereo_depth * 1000.0, config.stereo_rate, config.stereo_mix * 100.0);
    }

    fmt::print("\n🎵 СПЕКТР МОД (ТОЛЬКО ИЗ .dat):\n");
    fmt::print("   {:-<70}\n", "");
    fmt::print("   № | Частота(Гц) | Амплитуда | T60(с) | Атака(мс) | Название\n");
    fmt::print("   {:-<70}\n", "");
    for (size_t i = 0; i < modes.size() && i < 20; ++i) {
        fmt::print("   {:2d} | {:10.2f} | {:8.6f} | {:6.2f}  | {:7.1f}  | {}\n",
                   i+1, modes[i].frequency, modes[i].amplitude, 
                   modes[i].decay_time, modes[i].attack_time * 1000.0,
                   mode_names[i]);
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
    fmt::print("   Waveform: {}\n", 
               config.waveform_type == SynthesisConfig::WaveformType::SINE ? "SINE" :
               config.waveform_type == SynthesisConfig::WaveformType::TRIANGLE ? "TRIANGLE" :
               config.waveform_type == SynthesisConfig::WaveformType::SQUARE ? "SQUARE" :
               config.waveform_type == SynthesisConfig::WaveformType::SAWTOOTH ? "SAWTOOTH" : "MIXED");
    fmt::print("{:=^70}\n\n", "");
    
    return 0;
}