// =============================================================================
// synthesize.cpp - МОДАЛЬНЫЙ СИНТЕЗ КОЛОКОЛА
// =============================================================================
#include "bell_synthesis.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "ccx/ccx_parser.hpp"
#include "math_constants.hpp"
#include "core/physical_modal_selection.hpp"
#include "simple_json.hpp"



using namespace BellOptimizer;

int command_synthesize(int argc, char* argv[]) {
    std::cout << "\n======================================================================\n";
    std::cout << "======================================================================\n🔔 СИНТЕЗ КОЛОКОЛА (PHYSICAL MODAL SINE)\n";
    std::cout << "======================================================================\n\n";

    std::string filename = "best_bell.dat";
    bool filename_set = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: bell synth [file.dat] [--workdir <dir>]\n";
            std::cout << "Default input: best_bell.dat\n";
            return 0;
        }
        if (arg == "--workdir") {
            ++i;
            continue;
        }
        if (arg.rfind("--workdir=", 0) == 0) {
            continue;
        }
        if (!filename_set) {
            filename = arg;
            filename_set = true;
        } else {
            std::cerr << "Error: unknown synth option: " << arg << "\n";
            return 2;
        }
    }

    ModalData data = CCXParser::parse_file(filename);
    
    if (data.modes.empty()) {
        std::cerr << "❌ Ошибка: не найдены моды в файле " << filename << "\n";
        return 1;
    }

    // ============================================================
    // ЗАГРУЗКА КОНФИГУРАЦИИ ИЗ JSON
    // ============================================================
    SynthesisConfig config;
    
    try {
        std::ifstream config_file("synthesis.json");
        if (config_file.is_open()) {
            const auto cfg = bell::json::parse_file("synthesis.json");

            config.sample_rate = bell::json::integer(cfg, "sample_rate", 48000);
            config.oversample_factor = bell::json::integer(cfg, "oversample_factor", 4);
            config.duration = bell::json::number(cfg, "duration", 0.0);
            config.fade_in = bell::json::number(cfg, "fade_in", 0.003);
            config.fade_out = bell::json::number(cfg, "fade_out", 0.3);
            config.loss_factor = bell::json::number(cfg, "loss_factor", 0.0005);
            config.tail_factor = bell::json::number(cfg, "tail_factor", 1.2);
            config.target_peak = bell::json::number(cfg, "target_peak", 0.95);
            config.envelope_threshold = bell::json::number(cfg, "envelope_threshold", 0.00001);
            config.use_tail_trimming = bell::json::boolean(cfg, "use_tail_trimming", true);
            config.reverb_mix = bell::json::number(cfg, "reverb_mix", 0.0);
            config.reverb_feedback = bell::json::number(cfg, "reverb_feedback", 0.3);
            config.max_frequency = bell::json::number(cfg, "max_frequency", 20000.0);

            std::cout << "✅ Загружена конфигурация из synthesis.json\n";
        } else {
            std::cout << "⚠️ Файл synthesis.json не найден, используются значения по умолчанию\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Ошибка чтения synthesis.json: " << e.what() << "\n";
        std::cout << "   Используются значения по умолчанию\n";
    }

    // Определяем основную частоту
    double prime_frequency = 0.0;
    if (data.modes.size() > 1) {
        prime_frequency = data.modes[1].frequency_hz;
    } else {
        prime_frequency = data.modes[0].frequency_hz;
    }
    config.target_fundamental = prime_frequency;

    // Frequencies from the modal file are already absolute physical Hz.
    // Select by the requested synthesis band, not by a fixed mode count.
    PhysicalModalSelectionConfig selection_config;
    selection_config.maximum_frequency_hz = config.max_frequency;
    selection_config.sample_rate_hz = static_cast<double>(config.sample_rate);
    selection_config.nyquist_margin = 0.98;
    const auto selection = selectPhysicalSynthesisModes(data, selection_config);
    data.modes = selection.selected_modes;
    if (data.modes.empty() || prime_frequency > selection.effective_maximum_frequency_hz) {
        std::cerr << "No usable modes remain in the physical synthesis band, or the prime frequency is above its ceiling.\n";
        return 1;
    }
    std::cout << "   Physical modal band: 0-" << selection.effective_maximum_frequency_hz
              << " Hz; kept " << selection.selected_count << "/" << selection.input_count
              << " modes (outside band=" << selection.outside_frequency_band
              << ", numerical rejects=" << selection.numerical_rejections
              << ", weak retained=" << selection.weak_modes << ")\n";

    // ============================================================
    // ПОДГОТОВКА МОД
    // ============================================================
    std::vector<Mode> modes;
    std::vector<double> frequencies;
    std::vector<std::string> mode_names;
    
    // Acoustic transfer for the current generalized rotational excitation
    // followed by far-field monopole radiation. The excitation is a rigid-body
    // base rotation theta(t), with q_rotation=(-z,r). For an angular-velocity
    // impulse J_theta = integral(theta_ddot dt), the modal coordinate is
    //   q_m(t) = J_theta * part_rz / omega * sin(omega*t).
    // Therefore the signed pressure coefficient is
    //   C_m = -J_theta * part_rz * Qd * omega.
    // J_theta has units 1/s (angular velocity; radians are dimensionless).
    // It is NOT a force impulse or a torque impulse. The common J_theta, air
    // density and listener distance cancel during peak normalization. The
    // SIGN must be retained: part_rz and Qd both change sign when an
    // eigenvector is negated, so their product is invariant.
    double max_acoustic_excitation = 0.0;
    bool have_acoustic_coefficients = false;
    for (const auto& d : data.modes) {
        const double omega = 2.0 * bell::kPi * d.frequency_hz;
        const double pressure_coefficient =
            -d.part_rz * d.acoustic_volume_displacement * omega;
        if (std::abs(d.acoustic_volume_displacement) > 1e-14) {
            have_acoustic_coefficients = true;
            max_acoustic_excitation = std::max(
                max_acoustic_excitation, std::abs(pressure_coefficient));
        }
    }
    if (max_acoustic_excitation < 1e-14) {
        // Compatibility fallback for legacy CCX files without mode shapes.
        for (const auto& d : data.modes) {
            max_acoustic_excitation = std::max(max_acoustic_excitation,
                                               std::sqrt(std::max(0.0, d.effective_mass_norm)));
        }
    }
    if (max_acoustic_excitation < 1e-14) max_acoustic_excitation = 1.0;
    
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
        
        const double omega = 2.0 * bell::kPi * m.frequency;
        const double pressure_coefficient =
            -data.modes[i].part_rz * data.modes[i].acoustic_volume_displacement * omega;
        if (have_acoustic_coefficients) {
            const double magnitude = std::abs(pressure_coefficient);
            m.amplitude = magnitude / max_acoustic_excitation;
            // A negative signed pressure coefficient is exactly a pi phase
            // reversal relative to the common zero-phase impulse response.
            m.phase = (pressure_coefficient < 0.0) ? bell::kPi : 0.0;
        } else {
            m.amplitude = std::sqrt(std::max(0.0, data.modes[i].effective_mass_norm))
                        / max_acoustic_excitation;
            m.phase = 0.0;
        }
        
        if (m.frequency > 1.0) {
            // For structural loss factor eta, alpha ~= pi*f*eta and
            // T60 = ln(1000)/alpha for amplitude decay.
            m.decay_time = std::log(1000.0) / (bell::kPi * m.frequency * config.loss_factor);
        } else {
            m.decay_time = 5.0;
        }
        
        double freq_ratio_attack = m.frequency / config.target_fundamental;
        m.attack_time = 0.002 + 0.005 / std::sqrt(freq_ratio_attack);
        
        // The absolute generalized excitation impulse is intentionally absent from
        // this normalized synthesizer: a common J_theta scales every mode and
        // is removed by final peak normalization. It must not affect decay or
        // attack times. Those are determined by structural loss and the
        // synthesis envelope.

        m.attack_time = std::clamp(m.attack_time, 0.001, 0.03);
        m.decay_time = std::clamp(m.decay_time, 0.1, 30.0);
        
        modes.push_back(m);
    }

    if (modes.empty()) {
        std::cerr << "❌ Ошибка: нет мод для синтеза\n";
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

    std::cout << "\n📋 ПАРАМЕТРЫ:\n";
    std::cout << std::fixed << std::setprecision(1) << "   Основная частота (Prime): " << config.target_fundamental << " Гц\n";
    std::cout << "   Количество мод: " << modes.size() << "\n";
    std::cout << std::scientific << std::setprecision(2) << "   Коэффициент потерь: " << config.loss_factor << "\n" << std::defaultfloat;
    std::cout << std::fixed << std::setprecision(1) << "   Длительность: " << config.duration << " с\n";
    std::cout << std::fixed << std::setprecision(0) << "   Реверберация: " << config.reverb_mix * 100.0 << "%\n";

    std::cout << "\n🎵 СПЕКТР МОД (ТОЛЬКО ИЗ .dat):\n";
    std::cout << "   " << std::string(70, '-') << "\n";
    std::cout << "   № | Частота(Гц) | Амплитуда | T60(с) | Атака(мс) | Название\n";
    std::cout << "   " << std::string(70, '-') << "\n";
    for (size_t i = 0; i < modes.size() && i < 20; ++i) {
        std::cout << std::right << std::setw(5) << i + 1
                  << " | " << std::setw(10) << std::fixed << std::setprecision(2) << modes[i].frequency
                  << " | " << std::setw(8) << std::setprecision(6) << modes[i].amplitude
                  << " | " << std::setw(6) << std::setprecision(2) << modes[i].decay_time
                  << "  | " << std::setw(7) << std::setprecision(1) << modes[i].attack_time * 1000.0
                  << "  | " << mode_names[i] << "\n";
    }
    if (modes.size() > 20) {
        std::cout << "   ... и ещё " << modes.size() - 20 << " мод\n";
    }
    std::cout << "   " << std::string(70, '-') << "\n";

    auto start = std::chrono::high_resolution_clock::now();
    
    const long long frequency_tag = std::llround(config.target_fundamental);
    config.output_filename = "bell_" + std::to_string(frequency_tag) + "hz.wav";
    
    BellSynthesizer synth(config);
    synth.generate_modes(geo, modes);
    synth.synthesize();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\n   Время синтеза: " << ms.count() << " мс\n";
    synth.save_wav(config.output_filename);
    synth.analyze_spectrum(frequencies);
    synth.print_spectrum(modes);
    
    synth.save_all_envelopes_csv("envelopes_all.csv");
    
    std::cout << "\n   Выполняется FFT анализ...\n";
    
    size_t fft_size = synth.get_audio().size();
    synth.analyze_fft(fft_size);
    
    synth.print_fft_peaks(frequencies);
    synth.save_spectrum_csv("spectrum_" + std::to_string(frequency_tag) + "hz.csv", frequencies);
    
    std::cout << "\n======================================================================\n";
    std::cout << "✅ ГОТОВО: " << config.output_filename << "\n";
    std::cout << "   Источник: " << filename << "\n";
    std::cout << "   Мод: " << modes.size() << "\n";
    std::cout << std::fixed << std::setprecision(1) << "   Длительность: " << synth.get_actual_duration() << " с\n";
    std::cout << "   CSV спектр: spectrum_" << frequency_tag << "hz.csv\n";
    std::cout << "   CSV огибающих: envelopes_all.csv\n";
    std::cout << "======================================================================\n\n";
    
    return 0;
}