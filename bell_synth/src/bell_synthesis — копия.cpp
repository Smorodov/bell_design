// =============================================================================
// bell_synthesis.cpp - СИНТЕЗ ЗВУКА КОЛОКОЛА
// =============================================================================
#include "bell_synthesis.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iomanip>
#include <chrono>
#include <fstream>

#include <fmt/format.h>

#ifdef USE_FFTW
    #include <fftw3.h>
#endif

namespace BellOptimizer {

// =========================================================================
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
// =========================================================================

static double smoothstep(double t, double a, double b) {
    if (t <= a) return 0.0;
    if (t >= b) return 1.0;
    double x = (t - a) / (b - a);
    return x * x * (3.0 - 2.0 * x);
}

// =========================================================================
// GEOMETRYPARAMS
// =========================================================================

double GeometryParams::get_radius_at(double height_pos) const {
    if (!has_profile || profile_radius.empty()) {
        double t = std::max(0.0, std::min(1.0, height_pos));
        return r_min + (r_max - r_min) * t;
    }
    
    double pos = height_pos * (profile_radius.size() - 1);
    pos = std::clamp(pos, 0.0, static_cast<double>(profile_radius.size() - 1));
    size_t idx = static_cast<size_t>(pos);
    double frac = pos - idx;
    
    if (idx >= profile_radius.size() - 1) return profile_radius.back();
    return profile_radius[idx] * (1.0 - frac) + profile_radius[idx + 1] * frac;
}

double GeometryParams::get_thickness_at(double height_pos) const {
    if (!has_profile || profile_thickness.empty()) {
        return thickness;
    }
    
    double pos = height_pos * (profile_thickness.size() - 1);
    pos = std::clamp(pos, 0.0, static_cast<double>(profile_thickness.size() - 1));
    size_t idx = static_cast<size_t>(pos);
    double frac = pos - idx;
    
    if (idx >= profile_thickness.size() - 1) return profile_thickness.back();
    return profile_thickness[idx] * (1.0 - frac) + profile_thickness[idx + 1] * frac;
}

bool GeometryParams::is_monotonic() const {
    if (!has_profile || profile_radius.size() < 2) {
        return r_max >= r_min;
    }
    
    for (size_t i = 1; i < profile_radius.size(); ++i) {
        if (profile_radius[i] < profile_radius[i-1]) return false;
    }
    return true;
}

double GeometryParams::get_max_slope() const {
    if (!has_profile || profile_radius.size() < 2) {
        return (r_max - r_min) / height;
    }
    
    double max_slope = 0.0;
    double dh = height / (profile_radius.size() - 1);
    
    for (size_t i = 1; i < profile_radius.size(); ++i) {
        double slope = (profile_radius[i] - profile_radius[i-1]) / dh;
        if (std::abs(slope) > max_slope) max_slope = std::abs(slope);
    }
    return max_slope;
}

// =========================================================================
// КОНСТРУКТОР
// =========================================================================

BellSynthesizer::BellSynthesizer(const SynthesisConfig& config)
    : config_(config) {
    if (config_.duration > 0.0) {
        audio_.reserve(static_cast<size_t>(config_.duration * config_.sample_rate));
    }
}

// =========================================================================
// ВЫБОР ФОРМЫ ВОЛНЫ
// =========================================================================

PolyBLEP::Waveform BellSynthesizer::select_waveform(double frequency) const {
    switch (config_.waveform_type) {
        case SynthesisConfig::WaveformType::SINE:
            return PolyBLEP::SINE;
        case SynthesisConfig::WaveformType::TRIANGLE:
            return PolyBLEP::TRIANGLE;
        case SynthesisConfig::WaveformType::SQUARE:
            return PolyBLEP::SQUARE;
        case SynthesisConfig::WaveformType::SAWTOOTH:
            return PolyBLEP::SAWTOOTH;
        case SynthesisConfig::WaveformType::MIXED:
        default:
            if (frequency < 500.0) return PolyBLEP::TRIANGLE;
            else if (frequency < 2000.0) return PolyBLEP::SINE;
            else return PolyBLEP::SAWTOOTH;
    }
}

// =========================================================================
// РАСЧЕТ ОГИБАЮЩЕЙ
// =========================================================================

double BellSynthesizer::calculate_decay_envelope(double time, const Mode& mode) const {
    double t60 = std::max(mode.decay_time, 0.001);
    double attack_time = std::max(mode.attack_time, 0.001);
    
    // Атака
    double attack_end = attack_time * 1.5;
    double attack = std::min(1.0, time / attack_end);
    
    // Нелинейное демпфирование
    double amplitude = mode.amplitude * attack;
    double nonlinear_factor = 1.0 + config_.nonlinear_damping * amplitude * amplitude;
    double effective_t60 = t60 / nonlinear_factor;
    
    double freq_ratio = std::clamp(mode.frequency / 1000.0, 0.1, 10.0);
    
    // Три экспоненты с разными временами затухания
    double tau1 = effective_t60 / (1.0 + 0.5 * freq_ratio);
    double tau2 = effective_t60;
    double tau3 = effective_t60 * (1.0 + 0.3 / freq_ratio);
    
    tau1 = std::max(tau1, 0.001);
    if (tau1 > tau2) tau1 = tau2 * 0.5;
    if (tau3 < tau2) tau3 = tau2 * 1.5;
    tau3 = std::min(tau3, 30.0);
    
    // Точки переключения между экспонентами
    double t1 = effective_t60 * 0.05;
    double t2 = effective_t60 * 0.25;
    double t3 = effective_t60 * 0.50;
    double t4 = effective_t60 * 0.75;
    
    double e1 = 0.0, e2 = 0.0, e3 = 0.0;
    double w1 = 0.0, w2 = 0.0, w3 = 0.0;
    
    if (time < t1) {
        e1 = std::exp(-3.0 * time / tau1);
        w1 = 1.0;
    } else if (time < t2) {
        double alpha = smoothstep(time, t1, t2);
        w1 = 1.0 - alpha;
        w2 = alpha;
        e1 = std::exp(-3.0 * time / tau1);
        e2 = std::exp(-3.0 * time / tau2);
    } else if (time < t3) {
        w2 = 1.0;
        e2 = std::exp(-3.0 * time / tau2);
    } else if (time < t4) {
        double alpha = smoothstep(time, t3, t4);
        w2 = 1.0 - alpha;
        w3 = alpha;
        e2 = std::exp(-3.0 * time / tau2);
        double e2_at_t3 = std::exp(-3.0 * t3 / tau2);
        e3 = e2_at_t3 * std::exp(-3.0 * (time - t3) / tau3);
    } else {
        w3 = 1.0;
        double e2_at_t3 = std::exp(-3.0 * t3 / tau2);
        e3 = e2_at_t3 * std::exp(-3.0 * (time - t3) / tau3);
    }
    
    double envelope = w1 * e1 + w2 * e2 + w3 * e3;
    envelope *= attack;
    
    if (envelope < 0.0) envelope = 0.0;
    if (envelope < 1e-12) envelope = 0.0;
    
    return envelope;
}

// =========================================================================
// ГЕНЕРАЦИЯ МОДЫ С POLYBLEP
// =========================================================================

void BellSynthesizer::add_mode_with_polyblep(const Mode& mode, const EnvelopeData& envelope_data,
                                             std::vector<double>& buffer) {
    if (envelope_data.envelope.empty()) return;
    
    double dt = 1.0 / (config_.sample_rate * config_.oversample_factor);
    size_t max_samples = std::min(buffer.size(), 
                                  static_cast<size_t>(envelope_data.duration / dt));
    
    PolyBLEP osc(config_.sample_rate * config_.oversample_factor);
    osc.setFrequency(mode.frequency);
    osc.setWaveform(select_waveform(mode.frequency));
    osc.setAmplitude(1.0);
    
    double phase_norm = mode.phase / (2.0 * M_PI);
    phase_norm = phase_norm - std::floor(phase_norm);
    osc.sync(phase_norm);
    
    for (size_t i = 0; i < max_samples; ++i) {
        double t = i * dt;
        double envelope = get_envelope_value(envelope_data, t);
        double sample = osc.getAndInc();
        buffer[i] += mode.amplitude * envelope * sample;
    }
}

// =========================================================================
// ГЕНЕРАЦИЯ МОД
// =========================================================================

void BellSynthesizer::generate_modes(const GeometryParams& geo, 
                                     const std::vector<Mode>& base_modes) {
    geometry_ = geo;
    modes_ = base_modes;
    printf("   Modes for synthesis: %zu\n", modes_.size());
}

// =========================================================================
// РАСЧЕТ ЕСТЕСТВЕННОЙ ДЛИТЕЛЬНОСТИ
// =========================================================================

double BellSynthesizer::calculate_natural_duration() const {
    double max_t60 = 0.0;
    for (const auto& mode : modes_) {
        if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
    }
    
    max_t60 = std::max(max_t60, 0.001);
    double natural_duration = std::clamp(max_t60 * config_.tail_factor, 0.5, 30.0);
    
    return natural_duration;
}

// =========================================================================
// ПРЕДВАРИТЕЛЬНЫЙ РАСЧЕТ ОГИБАЮЩИХ (ОРИГИНАЛЬНАЯ РАБОЧАЯ ВЕРСИЯ)
// =========================================================================

void BellSynthesizer::precompute_envelopes(double duration) {
    envelope_cache_.clear();
    
    double table_step = std::max(config_.envelope_table_step, 0.0001);
    size_t table_size = std::max(static_cast<size_t>(duration / table_step) + 1, size_t(10));
    
    printf("   Computing envelopes (step %.1f ms, %zu points)...\n", 
           table_step * 1000.0, table_size);
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        EnvelopeData data;
        data.step = table_step;
        data.duration = duration;
        data.envelope.resize(table_size);
        
        for (size_t i = 0; i < table_size; ++i) {
            double t = i * table_step;
            if (t > duration) t = duration;
            data.envelope[i] = calculate_decay_envelope(t, modes_[m]);
        }
        
        envelope_cache_[m] = std::move(data);
    }
    
    printf("   Envelopes computed: %zu modes\n", modes_.size());
}

// =========================================================================
// ПОЛУЧЕНИЕ ЗНАЧЕНИЯ ОГИБАЮЩЕЙ (ОРИГИНАЛЬНАЯ РАБОЧАЯ ВЕРСИЯ)
// =========================================================================

double BellSynthesizer::get_envelope_value(const EnvelopeData& data, double time) const {
    if (data.envelope.empty()) return 0.0;
    if (time < 0.0) return data.envelope[0];
    
    if (time >= data.duration) {
        double t_extra = time - data.duration;
        double decay = std::exp(-3.0 * t_extra / data.duration);
        return data.envelope.back() * decay;
    }
    
    double idx_double = time / data.step;
    size_t idx = static_cast<size_t>(idx_double);
    double frac = idx_double - idx;
    
    if (idx >= data.envelope.size() - 1) return data.envelope.back();
    
    return data.envelope[idx] * (1.0 - frac) + data.envelope[idx + 1] * frac;
}

// =========================================================================
// ДАУНСЭМПЛИНГ
// =========================================================================

void BellSynthesizer::downsample_audio() {
    size_t total_samples = static_cast<size_t>(actual_duration_ * config_.sample_rate);
    if (total_samples < 1) total_samples = 1;
    
    audio_.assign(total_samples, 0.0);
    
    for (size_t i = 0; i < total_samples; ++i) {
        double sum = 0.0;
        size_t start = i * config_.oversample_factor;
        size_t end = std::min(start + config_.oversample_factor, audio_oversampled_.size());
        
        for (size_t j = start; j < end; ++j) {
            sum += audio_oversampled_[j];
        }
        audio_[i] = sum / (end - start);
    }
}

// =========================================================================
// ЭФФЕКТ: РЕВЕРБЕРАЦИЯ
// =========================================================================

void BellSynthesizer::add_reverb() {
    if (config_.reverb_mix < 0.001) return;
    
    size_t total_samples = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    double max_t60 = 0.0;
    for (const auto& mode : modes_) {
        if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
    }
    max_t60 = std::max(max_t60, 0.5);
    
    double reverb_rt60 = max_t60 * 0.6;
    double decay_per_sample = std::exp(-3.0 * dt / reverb_rt60);
    
    const int NUM_LINES = 3;
    std::vector<int> delays = {
        static_cast<int>(0.02 * config_.sample_rate),
        static_cast<int>(0.035 * config_.sample_rate),
        static_cast<int>(0.055 * config_.sample_rate)
    };
    
    std::vector<std::vector<double>> delay_lines(NUM_LINES);
    std::vector<size_t> positions(NUM_LINES, 0);
    
    for (int i = 0; i < NUM_LINES; ++i) {
        delay_lines[i].assign(delays[i] + 1, 0.0);
    }
    
    std::vector<double> reverb(total_samples, 0.0);
    double mix = std::clamp(config_.reverb_mix, 0.0, 0.5);
    double feedback = std::clamp(config_.reverb_feedback, 0.0, 0.9);
    
    for (size_t i = 0; i < total_samples; ++i) {
        double wet = 0.0;
        double input = audio_[i] * 0.5;
        
        for (int ch = 0; ch < NUM_LINES; ++ch) {
            size_t pos = positions[ch];
            double delayed = delay_lines[ch][pos];
            double output = delayed * decay_per_sample;
            wet += output * 0.33;
            delay_lines[ch][pos] = input + output * feedback * 0.7;
            positions[ch] = (pos + 1) % delay_lines[ch].size();
        }
        
        reverb[i] = wet;
    }
    
    for (size_t i = 0; i < total_samples; ++i) {
        audio_[i] = audio_[i] * (1.0 - mix) + reverb[i] * mix;
    }
}

// =========================================================================
// ЭФФЕКТ: FADE
// =========================================================================

void BellSynthesizer::apply_fade(double fade_in_time, double fade_out_time) {
    size_t total_samples = audio_.size();
    
    size_t fade_in_samples = static_cast<size_t>(fade_in_time * config_.sample_rate);
    size_t fade_out_samples = static_cast<size_t>(fade_out_time * config_.sample_rate);
    fade_out_samples = std::min(fade_out_samples, total_samples / 4);
    
    for (size_t i = 0; i < fade_in_samples && i < total_samples; ++i) {
        double t = static_cast<double>(i) / fade_in_samples;
        double fade = t * t * (3.0 - 2.0 * t);
        audio_[i] *= fade;
    }
    
    for (size_t i = 0; i < fade_out_samples && i < total_samples; ++i) {
        size_t idx = total_samples - 1 - i;
        double t = static_cast<double>(i) / fade_out_samples;
        double fade = 1.0 - t * t * (3.0 - 2.0 * t);
        audio_[idx] *= fade;
    }
}

// =========================================================================
// НОРМАЛИЗАЦИЯ
// =========================================================================

void BellSynthesizer::normalize(double peak_level) {
    if (audio_.empty()) return;
    
    peak_amplitude_ = 0.0;
    double rms_sum = 0.0;
    
    for (double sample : audio_) {
        double abs_val = std::abs(sample);
        if (abs_val > peak_amplitude_) peak_amplitude_ = abs_val;
        rms_sum += sample * sample;
    }
    
    rms_amplitude_ = std::sqrt(rms_sum / audio_.size());
    
    if (peak_amplitude_ < 1e-10) return;
    
    double gain = peak_level / peak_amplitude_;
    if (config_.use_rms_normalization && rms_amplitude_ > 1e-10) {
        double rms_gain = config_.target_rms / rms_amplitude_;
        gain = std::min(gain, rms_gain * 1.2);
    }
    
    for (double& sample : audio_) {
        sample *= gain;
    }
}

// =========================================================================
// ОБРЕЗКА ТИШИНЫ
// =========================================================================

void BellSynthesizer::trim_silence(double threshold) {
    if (audio_.empty()) return;
    if (!config_.use_tail_trimming) return;
    
    // Проверяем начало
    size_t start_sample = 0;
    for (size_t i = 0; i < audio_.size(); ++i) {
        if (std::abs(audio_[i]) > threshold) {
            start_sample = i;
            break;
        }
    }
    
    if (start_sample == 0) {
        double max_amp = 0.0;
        for (size_t i = 0; i < std::min(size_t(1000), audio_.size()); ++i) {
            max_amp = std::max(max_amp, std::abs(audio_[i]));
        }
        if (max_amp < threshold) {
            printf("   Skipping trim: entire signal is quiet (max_amp=%.6f < threshold=%.6f)\n", 
                   max_amp, threshold);
            return;
        }
        start_sample = 0;
    }
    
    if (start_sample > 0) {
        size_t new_start = start_sample - 1;
        if (new_start > 0) {
            std::copy(audio_.begin() + new_start, audio_.end(), audio_.begin());
            audio_.resize(audio_.size() - new_start);
            printf("   Trimmed silence at start: %zu samples\n", new_start);
        }
    }
    
    // Проверяем конец
    bool has_silence_at_end = true;
    size_t check_samples = std::min(size_t(1000), audio_.size());
    for (size_t i = audio_.size() - check_samples; i < audio_.size(); ++i) {
        if (std::abs(audio_[i]) > threshold) {
            has_silence_at_end = false;
            break;
        }
    }
    
    if (!has_silence_at_end) {
        printf("   No silence at end, skipping tail trim\n");
        return;
    }
    
    // Ищем 50 мс тишины
    double silence_duration = 0.05;
    size_t silence_samples = static_cast<size_t>(silence_duration * config_.sample_rate);
    size_t min_silence_samples = std::max(silence_samples, size_t(10));
    
    size_t silence_start = 0;
    size_t silence_count = 0;
    bool found = false;
    
    for (size_t i = audio_.size(); i > 0; --i) {
        size_t idx = i - 1;
        if (std::abs(audio_[idx]) < threshold) {
            silence_count++;
            if (silence_count >= min_silence_samples) {
                silence_start = idx;
                found = true;
                break;
            }
        } else {
            silence_count = 0;
        }
    }
    
    if (!found) {
        printf("   No silence > %.0f ms found, skipping trim\n", silence_duration * 1000.0);
        return;
    }
    
    size_t new_size = silence_start + 1;
    size_t min_size = static_cast<size_t>(0.5 * config_.sample_rate);
    if (new_size < min_size) {
        printf("   Skip trim: would leave only %zu samples (< 0.5s)\n", new_size);
        return;
    }
    
    size_t padding = static_cast<size_t>(0.01 * config_.sample_rate);
    if (new_size > padding) {
        new_size -= padding;
    }
    
    if (new_size < audio_.size() && new_size > 100) {
        double old_duration = static_cast<double>(audio_.size()) / config_.sample_rate;
        double new_duration = static_cast<double>(new_size) / config_.sample_rate;
        
        printf("   Trimming tail: %.2f s -> %.2f s (silence > %.0f ms, threshold=%.6f)\n",
               old_duration, new_duration, silence_duration * 1000.0, threshold);
        
        size_t fade_samples = static_cast<size_t>(0.02 * config_.sample_rate);
        fade_samples = std::min(fade_samples, new_size / 4);
        
        if (fade_samples > 0) {
            for (size_t i = 0; i < fade_samples && i < new_size; ++i) {
                double t = static_cast<double>(i) / fade_samples;
                double fade = 0.5 * (1.0 + std::cos(M_PI * (1.0 - t)));
                size_t idx = new_size - 1 - i;
                audio_[idx] *= fade;
            }
        }
        
        audio_.resize(new_size);
    }
}

// =========================================================================
// ОСНОВНОЙ СИНТЕЗ
// =========================================================================

void BellSynthesizer::synthesize() {
    if (modes_.empty()) {
        fmt::print(stderr, "Error: no modes for synthesis\n");
        return;
    }
    
    double natural_duration = calculate_natural_duration();
    printf("\n=== DURATION CALC ===\n");
    printf("   Natural (by T60): %.2f s\n", natural_duration);
    
    if (config_.duration <= 0.0) {
        actual_duration_ = natural_duration;
        printf("   Mode: natural duration\n");
    } else {
        actual_duration_ = std::min(natural_duration, std::clamp(config_.duration, 0.1, 30.0));
        printf("   Requested: %.2f s\n", config_.duration);
        printf("   Final: %.2f s\n", actual_duration_);
    }
    
    printf("   actual_duration_ = %.6f s\n", actual_duration_);
    
    precompute_envelopes(actual_duration_);
    
    size_t total_samples_os = static_cast<size_t>(actual_duration_ * 
                                                   config_.sample_rate * 
                                                   config_.oversample_factor);
    printf("   total_samples_os = %zu\n", total_samples_os);
    
    if (total_samples_os < 1) total_samples_os = 1;
    
    audio_oversampled_.assign(total_samples_os, 0.0);
    
    printf("\n   Synthesizing signal (PolyBLEP, %zu modes)...\n", modes_.size());
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        auto it = envelope_cache_.find(m);
        if (it != envelope_cache_.end()) {
            add_mode_with_polyblep(modes_[m], it->second, audio_oversampled_);
        }
    }
    
    // Диагностика: проверяем audio_oversampled_
    double max_os = 0.0;
    size_t check_os = std::min(size_t(1000), audio_oversampled_.size());
    for (size_t i = 0; i < check_os; ++i) {
        max_os = std::max(max_os, std::abs(audio_oversampled_[i]));
    }
    printf("   audio_oversampled_ size = %zu, max amp (first 1000) = %.6f\n", 
           audio_oversampled_.size(), max_os);
    
    downsample_audio();
    
    // Диагностика: проверяем audio_ после даунсэмплинга
    double max_audio = 0.0;
    size_t check_audio = std::min(size_t(1000), audio_.size());
    for (size_t i = 0; i < check_audio; ++i) {
        max_audio = std::max(max_audio, std::abs(audio_[i]));
    }
    printf("   audio_ size = %zu, max amp (first 1000) = %.6f\n", 
           audio_.size(), max_audio);
    
    double signal_duration = static_cast<double>(audio_.size()) / static_cast<double>(config_.sample_rate);
    printf("   Signal: %zu samples, %.2f s\n", audio_.size(), signal_duration);
    
    // ============================================================
    // ПРИМЕНЕНИЕ ЭФФЕКТОВ (ТОЛЬКО РЕВЕРБЕРАЦИЯ)
    // ============================================================
    
    printf("\n   Applying effects...\n");
    
    if (config_.reverb_mix > 0.001) {
        printf("     - Reverb (mix=%.2f)\n", config_.reverb_mix);
        add_reverb();
    }
    
    // Диагностика: после реверберации
    double max_reverb = 0.0;
    size_t check_reverb = std::min(size_t(1000), audio_.size());
    for (size_t i = 0; i < check_reverb; ++i) {
        max_reverb = std::max(max_reverb, std::abs(audio_[i]));
    }
    printf("   audio_ after reverb size = %zu, max amp = %.6f\n", 
           audio_.size(), max_reverb);
    
    // ============================================================
    // FADE И ОБРЕЗКА
    // ============================================================
    
    apply_fade(config_.fade_in, config_.fade_out);
    
    // Диагностика: после fade
    double max_fade = 0.0;
    size_t check_fade = std::min(size_t(1000), audio_.size());
    for (size_t i = 0; i < check_fade; ++i) {
        max_fade = std::max(max_fade, std::abs(audio_[i]));
    }
    printf("   audio_ after fade size = %zu, max amp = %.6f\n", 
           audio_.size(), max_fade);
    
    if (config_.envelope_threshold > 0.0 && config_.use_tail_trimming) {
        printf("   Trimming silence (threshold=%.6f)...\n", config_.envelope_threshold);
        trim_silence(config_.envelope_threshold);
    }
    
    // Диагностика: после обрезки
    double max_trim = 0.0;
    size_t check_trim = std::min(size_t(1000), audio_.size());
    for (size_t i = 0; i < check_trim; ++i) {
        max_trim = std::max(max_trim, std::abs(audio_[i]));
    }
    printf("   audio_ after trim size = %zu, max amp = %.6f\n", 
           audio_.size(), max_trim);
    
    normalize(config_.target_peak);
    
    double final_duration = static_cast<double>(audio_.size()) / static_cast<double>(config_.sample_rate);
    printf("\n✅ Synthesis complete\n");
    printf("   Final duration: %.2f s\n", final_duration);
}

// =========================================================================
// СОХРАНЕНИЕ WAV
// =========================================================================

bool BellSynthesizer::save_wav(const std::string& filename) {
    if (audio_.empty()) {
        fmt::print(stderr, "Error: no audio to save\n");
        return false;
    }
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        fmt::print(stderr, "Error: failed to create file {}\n", filename);
        return false;
    }
    
    std::vector<int16_t> pcm_audio(audio_.size());
    for (size_t i = 0; i < audio_.size(); ++i) {
        double sample = std::max(-1.0, std::min(1.0, audio_[i]));
        pcm_audio[i] = static_cast<int16_t>(sample * 32767.0);
    }
    
    int data_size = static_cast<int>(pcm_audio.size() * sizeof(int16_t));
    write_wav_header(file, data_size);
    file.write(reinterpret_cast<const char*>(pcm_audio.data()), data_size);
    
    if (!file.good()) {
        fmt::print(stderr, "Error: write failed\n");
        file.close();
        return false;
    }
    
    file.close();
    printf("   Saved: %s\n", filename.c_str());
    printf("   Size: %d KB\n", data_size / 1024);
    
    return true;
}

void BellSynthesizer::write_wav_header(std::ofstream& file, int data_size) const {
    int sample_rate = config_.sample_rate;
    int channels = 1;
    int bits_per_sample = 16;
    
    int byte_rate = sample_rate * channels * bits_per_sample / 8;
    int block_align = channels * bits_per_sample / 8;
    int subchunk2_size = data_size;
    int chunk_size = 36 + subchunk2_size;
    
    file.write("RIFF", 4);
    file.write(reinterpret_cast<const char*>(&chunk_size), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    
    int subchunk1_size = 16;
    short audio_format = 1;
    short num_channels = static_cast<short>(channels);
    int sample_rate_int = sample_rate;
    short bits_per_sample_short = static_cast<short>(bits_per_sample);
    
    file.write(reinterpret_cast<const char*>(&subchunk1_size), 4);
    file.write(reinterpret_cast<const char*>(&audio_format), 2);
    file.write(reinterpret_cast<const char*>(&num_channels), 2);
    file.write(reinterpret_cast<const char*>(&sample_rate_int), 4);
    file.write(reinterpret_cast<const char*>(&byte_rate), 4);
    file.write(reinterpret_cast<const char*>(&block_align), 2);
    file.write(reinterpret_cast<const char*>(&bits_per_sample_short), 2);
    
    file.write("data", 4);
    file.write(reinterpret_cast<const char*>(&subchunk2_size), 4);
}

// =========================================================================
// АНАЛИЗ СПЕКТРА
// =========================================================================

void BellSynthesizer::analyze_spectrum(const std::vector<double>& frequencies) {
    if (audio_.empty()) return;
    
    size_t n = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    printf("\n--- SPECTRUM ---\n");
    printf("   Target frequencies (Hz):\n");
    printf("   ------------------------\n");
    
    double time_norm = 1.0 / n;
    
    for (double freq : frequencies) {
        if (freq > config_.max_frequency) continue;
        
        double real = 0.0, imag = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double t = i * dt;
            double phase = 2.0 * M_PI * freq * t;
            real += audio_[i] * std::cos(phase);
            imag += audio_[i] * std::sin(phase);
        }
        double amp = std::sqrt(real * real + imag * imag) * time_norm;
        
        if (amp > 0.0001) {
            printf("   %.1f Hz : %.3f mV\n", freq, amp * 1000.0);
        }
    }
    
    printf("   ------------------------\n");
}

void BellSynthesizer::print_spectrum(const std::vector<Mode>& modes) const {
    if (modes.empty()) return;
    
    printf("\n--- MODE SPECTRUM ---\n");
    printf("   # | Freq(Hz) | T60(s) | Attack(s) | Amplitude\n");
    printf("   ---------------------------------------------\n");
    
    for (size_t i = 0; i < modes.size() && i < 25; ++i) {
        const auto& m = modes[i];
        printf("   %2zu | %10.1f | %6.1f | %7.3f | %8.6f\n",
               i+1, m.frequency, m.decay_time, m.attack_time, m.amplitude);
    }
    printf("   ---------------------------------------------\n");
}

// =========================================================================
// FFT АНАЛИЗ
// =========================================================================

#ifdef USE_FFTW
void BellSynthesizer::analyze_fft(int fft_size) {
    if (audio_.empty()) {
        fmt::print(stderr, "Error: no audio for FFT analysis\n");
        return;
    }

    size_t max_fft_size = 65536;
    size_t target = (fft_size <= 0) ? std::min(audio_.size(), max_fft_size) 
                                    : std::min(static_cast<size_t>(fft_size), 
                                               std::min(audio_.size(), max_fft_size));
    
    size_t n = 1;
    while (n < target) n <<= 1;

    if (n < 2) {
        fmt::print(stderr, "Error: FFT size too small (n < 2)\n");
        return;
    }

    double resolution = static_cast<double>(config_.sample_rate) / static_cast<double>(n);
    printf("   Performing FFT (size: %zu, resolution: %.1f Hz)...\n", n, resolution);

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    
    if (!in || !out) {
        fmt::print(stderr, "Error: failed to allocate memory for FFT\n");
        if (in) fftw_free(in);
        if (out) fftw_free(out);
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        double window = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (n - 1)));
        if (i < target) {
            in[i][0] = audio_[i] * window;
            in[i][1] = 0.0;
        } else {
            in[i][0] = 0.0;
            in[i][1] = 0.0;
        }
    }

    fftw_plan plan = fftw_plan_dft_1d((int)n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    double dt = 1.0 / config_.sample_rate;
    fft_frequencies_.clear();
    fft_magnitudes_.clear();
    
    double window_compensation = 2.0;
    
    for (size_t i = 0; i < n / 2; ++i) {
        double freq = static_cast<double>(i) / (static_cast<double>(n) * dt);
        double mag = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / static_cast<double>(n) * window_compensation;
        fft_frequencies_.push_back(freq);
        fft_magnitudes_.push_back(mag);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}
#endif

void BellSynthesizer::print_fft_peaks(const std::vector<double>& target_frequencies) {
    if (fft_frequencies_.empty()) {
        fmt::print(stderr, "Please run analyze_fft() first\n");
        return;
    }

    printf("\n--- SPECTRUM PEAKS ---\n");
    printf("   Frequency(Hz) | Amplitude\n");
    printf("   -------------------------\n");

    double freq_resolution = (fft_frequencies_.size() > 1) ? 
                             (fft_frequencies_[1] - fft_frequencies_[0]) : 0.5;
    double tolerance = std::max(0.5, freq_resolution * 3.0);
    
    for (double target_freq : target_frequencies) {
        if (target_freq < 1.0) continue;
        
        double best_freq = 0.0;
        double best_mag = 0.0;
        
        for (size_t i = 0; i < fft_frequencies_.size(); ++i) {
            double diff = std::abs(fft_frequencies_[i] - target_freq);
            if (diff < tolerance && fft_magnitudes_[i] > best_mag) {
                best_freq = fft_frequencies_[i];
                best_mag = fft_magnitudes_[i];
            }
        }
        
        if (best_mag > 0.0) {
            printf("   %.2f Hz -> %.6f (found %.2f Hz)\n", 
                   target_freq, best_mag, best_freq);
        } else {
            printf("   %.2f Hz -> NOT FOUND\n", target_freq);
        }
    }
    
    printf("   -------------------------\n");
}

void BellSynthesizer::save_spectrum_csv(const std::string& filename, 
                                        const std::vector<double>& target_frequencies) {
    if (fft_frequencies_.empty()) {
        fmt::print(stderr, "Please run analyze_fft() first\n");
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        fmt::print(stderr, "Error: failed to create {}\n", filename);
        return;
    }

    file << "TargetFrequency(Hz),FoundFrequency(Hz),Magnitude\n";
    
    double freq_resolution = (fft_frequencies_.size() > 1) ? 
                             (fft_frequencies_[1] - fft_frequencies_[0]) : 0.5;
    double tolerance = std::max(0.5, freq_resolution * 3.0);
    int found_count = 0;
    
    for (double target_freq : target_frequencies) {
        if (target_freq < 1.0) continue;
        
        double best_freq = 0.0;
        double best_mag = 0.0;
        
        for (size_t i = 0; i < fft_frequencies_.size(); ++i) {
            double diff = std::abs(fft_frequencies_[i] - target_freq);
            if (diff < tolerance && fft_magnitudes_[i] > best_mag) {
                best_freq = fft_frequencies_[i];
                best_mag = fft_magnitudes_[i];
            }
        }
        
        if (best_mag > 0.0) {
            file << std::fixed << std::setprecision(2) 
                 << target_freq << ","
                 << std::setprecision(2) << best_freq << ","
                 << std::scientific << std::setprecision(6) << best_mag << "\n";
            found_count++;
        } else {
            file << std::fixed << std::setprecision(2) 
                 << target_freq << ",NOT_FOUND,0.0\n";
        }
    }
    
    file.close();
    
    printf("   Spectrum saved: %s (%d of %zu frequencies found)\n", 
           filename.c_str(), found_count, target_frequencies.size());
}

void BellSynthesizer::save_all_envelopes_csv(const std::string& filename) {
    if (envelope_cache_.empty()) {
        fmt::print(stderr, "Error: envelopes not computed yet\n");
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        fmt::print(stderr, "Error: failed to create {}\n", filename);
        return;
    }

    double dt = 1.0 / (config_.sample_rate * config_.oversample_factor);
    size_t total_samples = static_cast<size_t>(actual_duration_ * 
                                                config_.sample_rate * 
                                                config_.oversample_factor);
    double step = 0.01;
    int step_samples = static_cast<int>(step / dt);
    if (step_samples < 1) step_samples = 1;

    file << "Time(ms)";
    for (size_t i = 0; i < modes_.size(); ++i) {
        file << ",Mode" << (i+1) << "_" << std::fixed << std::setprecision(1) << modes_[i].frequency << "Hz";
    }
    file << "\n";

    for (int i = 0; i < (int)total_samples; i += step_samples) {
        double t = i * dt;
        file << std::fixed << std::setprecision(3) << (t * 1000.0);
        
        for (size_t m = 0; m < modes_.size(); ++m) {
            auto it = envelope_cache_.find(m);
            if (it != envelope_cache_.end()) {
                double envelope = get_envelope_value(it->second, t);
                file << "," << std::setprecision(6) << envelope;
            } else {
                file << ",0.0";
            }
        }
        file << "\n";
    }
    
    file.close();
    printf("   All envelopes saved: %s\n", filename.c_str());
}

} // namespace BellOptimizer