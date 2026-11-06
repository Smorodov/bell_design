// =============================================================================
// bell_synthesis.cpp - СИНТЕЗ ЗВУКА КОЛОКОЛА
// =============================================================================
#include "bell_synthesis.hpp"
#include "../audio/wav_writer.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iomanip>
#include <chrono>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "math_constants.hpp"



namespace BellOptimizer {

// =========================================================================
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
// =========================================================================

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
    if (!std::isfinite(config_.duration) || config_.duration <= 0.0 ||
        config_.sample_rate <= 0 || config_.oversample_factor <= 0) {
        throw std::invalid_argument("SynthesisConfig requires finite positive duration, sample rate, and oversample factor");
    }
    const long double base_samples = static_cast<long double>(config_.duration) * config_.sample_rate;
    // synthesize() clamps its actual generation window to [0.5, 30] seconds,
    // so validate against that maximum, not just the requested duration.
    const long double max_base_samples = 30.0L * config_.sample_rate;
    const long double max_oversampled_samples = max_base_samples * config_.oversample_factor;
    const auto max_samples = static_cast<long double>(std::vector<double>().max_size());
    const auto max_wav_samples = static_cast<long double>((std::numeric_limits<int>::max() - 36) / 2);
    if (!std::isfinite(base_samples) || !std::isfinite(max_base_samples) ||
        !std::isfinite(max_oversampled_samples) || base_samples > max_samples ||
        max_base_samples > max_wav_samples || max_oversampled_samples > max_samples) {
        throw std::length_error("SynthesisConfig sample count exceeds supported WAV/vector capacity");
    }
    audio_.reserve(static_cast<size_t>(base_samples));
}

// =========================================================================
// РАСЧЕТ ОГИБАЮЩЕЙ
// =========================================================================

double BellSynthesizer::calculate_decay_envelope(double time, const Mode& mode) const {
    const double t60 = std::max(mode.decay_time, 0.001);
    const double attack_time = std::max(mode.attack_time, 0.001);

    // T60 is an amplitude decay time: A(T60) = A0 / 1000.
    // Therefore A(t) = exp(-ln(1000) * t / T60).
    const double decay = std::exp(-std::log(1000.0) * time / t60);

    double attack = 1.0;
    if (time < attack_time) {
        const double x = std::clamp(time / attack_time, 0.0, 1.0);
        // Smooth attack without introducing an additional decay mechanism.
        attack = x * x * (3.0 - 2.0 * x);
    }

    const double envelope = attack * decay;
    return (envelope < 1e-12) ? 0.0 : envelope;
}

// =========================================================================
// ГЕНЕРАЦИЯ СИНУСОИДАЛЬНОЙ МОДЫ
// =========================================================================

void BellSynthesizer::add_mode(const Mode& mode, const EnvelopeData& envelope_data,
                                std::vector<double>& buffer) {
    if (envelope_data.envelope.empty()) return;
    
    double dt = 1.0 / (config_.sample_rate * config_.oversample_factor);
    size_t max_samples = std::min(buffer.size(), 
                                  static_cast<size_t>(envelope_data.duration / dt));
    
    double phase = mode.phase;
    const double phase_step = 2.0 * bell::kPi * mode.frequency * dt;

    for (size_t i = 0; i < max_samples; ++i) {
        const double t = i * dt;
        const double envelope = get_envelope_value(envelope_data, t);
        buffer[i] += mode.amplitude * envelope * std::sin(phase);
        phase += phase_step;
        if (phase > 2.0 * bell::kPi) {
            phase = std::fmod(phase, 2.0 * bell::kPi);
        }
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
// ПРЕДВАРИТЕЛЬНЫЙ РАСЧЕТ ОГИБАЮЩИХ (С ИНДИВИДУАЛЬНОЙ ДЛИТЕЛЬНОСТЬЮ)
// =========================================================================

void BellSynthesizer::precompute_envelopes(double global_duration) {
    envelope_cache_.clear();
    
    double table_step = std::max(config_.envelope_table_step, 0.0001);
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        // Индивидуальная длительность для каждой моды
        double mode_duration = modes_[m].decay_time * 3.0;
        mode_duration = std::clamp(mode_duration, 0.5, global_duration);
        
        size_t table_size = std::max(static_cast<size_t>(mode_duration / table_step) + 1, size_t(10));
        
        EnvelopeData data;
        data.step = table_step;
        data.duration = mode_duration;
        data.envelope.resize(table_size);
        
        for (size_t i = 0; i < table_size; ++i) {
            double t = i * table_step;
            if (t > mode_duration) t = mode_duration;
            data.envelope[i] = calculate_decay_envelope(t, modes_[m]);
        }
        
        // Плавное затухание к нулю на последних 10%
        size_t fade_start = static_cast<size_t>(table_size * 0.9);
        if (fade_start < table_size) {
            for (size_t i = fade_start; i < table_size; ++i) {
                double t = static_cast<double>(i - fade_start) / (table_size - fade_start);
                double fade = 1.0 - t * t * (3.0 - 2.0 * t);
                data.envelope[i] *= fade;
            }
        }
        
        envelope_cache_[m] = std::move(data);
    }
    
    printf("   Envelopes computed: %zu modes\n", modes_.size());
}

// =========================================================================
// РАСЧЕТ ФАКТИЧЕСКОЙ ДЛИТЕЛЬНОСТИ ПО СУММАРНОЙ ОГИБАЮЩЕЙ
// =========================================================================

double BellSynthesizer::calculate_actual_duration_from_envelopes(double threshold) {
    if (envelope_cache_.empty()) return 0.5;
    
    // Находим максимальную длительность среди всех огибающих
    double max_duration = 0.0;
    for (const auto& pair : envelope_cache_) {
        if (pair.second.duration > max_duration) {
            max_duration = pair.second.duration;
        }
    }
    
    if (max_duration < 0.01) return 0.5;
    
    // Шаг для сканирования (1 мс)
    double step = 0.001;
    size_t steps = static_cast<size_t>(max_duration / step) + 1;
    
    // Идём с конца и ищем момент, когда суммарная огибающая > порога
    for (size_t i = steps; i > 0; --i) {
        double time = i * step;
        if (time > max_duration) continue;
        
        double audible_envelope = 0.0;
        for (size_t m = 0; m < modes_.size(); ++m) {
            const auto it = envelope_cache_.find(m);
            if (it == envelope_cache_.end()) continue;
            audible_envelope = std::max(
                audible_envelope,
                std::abs(modes_[m].amplitude) * get_envelope_value(it->second, time));
        }

        // The audible tail is governed by the strongest remaining modal
        // contribution, not by an unweighted sum of all envelopes.
        if (audible_envelope > threshold) {
            // Добавляем небольшой запас (0.1 секунды)
            double result = time + 0.1;
            return std::min(result, max_duration);
        }
    }
    
    return max_duration;
}

// =========================================================================
// ПОЛУЧЕНИЕ ЗНАЧЕНИЯ ОГИБАЮЩЕЙ
// =========================================================================

double BellSynthesizer::get_envelope_value(const EnvelopeData& data, double time) const {
    if (data.envelope.empty()) return 0.0;
    if (time < 0.0) return data.envelope[0];
    
    // Выход за пределы длительности - возвращаем 0
    if (time >= data.duration) {
        return 0.0;
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
    const int oversample = std::max(config_.oversample_factor, 1);
    size_t total_samples = static_cast<size_t>(actual_duration_ * config_.sample_rate);
    if (total_samples < 1) total_samples = 1;

    audio_.assign(total_samples, 0.0);

    for (size_t i = 0; i < total_samples; ++i) {
        const size_t start = i * static_cast<size_t>(oversample);
        const size_t end = std::min(start + static_cast<size_t>(oversample),
                                    audio_oversampled_.size());
        if (start >= end) break;

        double sum = 0.0;
        for (size_t j = start; j < end; ++j) sum += audio_oversampled_[j];
        audio_[i] = sum / static_cast<double>(end - start);
    }
}

// =========================================================================
// ЭФФЕКТ: РЕВЕРБЕРАЦИЯ
// =========================================================================

void BellSynthesizer::add_reverb() {
    if (config_.reverb_mix < 0.001) return;
    
    size_t total_samples = audio_.size();
    
    double max_t60 = 0.0;
    for (const auto& mode : modes_) {
        if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
    }
    max_t60 = std::max(max_t60, 0.5);
    
    double reverb_rt60 = max_t60 * 0.6;
    const double dt = 1.0 / config_.sample_rate;
    double decay_per_sample = std::exp(-std::log(1000.0) * dt / reverb_rt60);
    
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
    const size_t total_samples = audio_.size();
    if (total_samples == 0) return;

    double fade_in_seconds = std::max(fade_in_time, 0.0);
    double fade_out_seconds = std::max(fade_out_time, 0.0);
    const double total_seconds = static_cast<double>(total_samples) / config_.sample_rate;
    const double requested = fade_in_seconds + fade_out_seconds;
    if (requested > total_seconds && requested > 0.0) {
        const double scale = total_seconds / requested;
        fade_in_seconds *= scale;
        fade_out_seconds *= scale;
    }

    const size_t fade_in_samples = std::min(
        total_samples, static_cast<size_t>(fade_in_seconds * config_.sample_rate));
    const size_t fade_out_samples = std::min(
        total_samples, static_cast<size_t>(fade_out_seconds * config_.sample_rate));

    for (size_t i = 0; i < fade_in_samples; ++i) {
        const double t = static_cast<double>(i) / fade_in_samples;
        audio_[i] *= t * t * (3.0 - 2.0 * t);
    }

    for (size_t i = 0; i < fade_out_samples; ++i) {
        const size_t idx = total_samples - 1 - i;
        const double t = static_cast<double>(i) / fade_out_samples;
        audio_[idx] *= 1.0 - t * t * (3.0 - 2.0 * t);
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
        peak_amplitude_ = std::max(peak_amplitude_, std::abs(sample));
        rms_sum += sample * sample;
    }
    rms_amplitude_ = std::sqrt(rms_sum / audio_.size());

    const double target = std::clamp(peak_level, 0.0, 1.0);
    if (peak_amplitude_ < 1e-12 || target <= 0.0) return;

    const double gain = target / peak_amplitude_;
    for (double& sample : audio_) sample *= gain;

    peak_amplitude_ = 0.0;
    rms_sum = 0.0;
    for (double sample : audio_) {
        peak_amplitude_ = std::max(peak_amplitude_, std::abs(sample));
        rms_sum += sample * sample;
    }
    rms_amplitude_ = std::sqrt(rms_sum / audio_.size());
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
                double fade = 0.5 * (1.0 + std::cos(bell::kPi * (1.0 - t)));
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
        std::cerr << "Error: no modes for synthesis\n";
        return;
    }
    
    double natural_duration = calculate_natural_duration();
    printf("\n=== DURATION CALC ===\n");
    printf("   Natural (by T60): %.2f s\n", natural_duration);
    
    double requested_duration = (config_.duration <= 0.0) ? natural_duration : config_.duration;
    double global_duration = std::clamp(requested_duration, 0.5, 30.0);
    printf("   Initial duration: %.2f s\n", global_duration);
    
    // ============================================================
    // 1. РАССЧИТЫВАЕМ ВСЕ ОГИБАЮЩИЕ
    // ============================================================
    
    precompute_envelopes(global_duration);
    
    // ============================================================
    // 2. ОПРЕДЕЛЯЕМ ФАКТИЧЕСКУЮ ДЛИТЕЛЬНОСТЬ ПО СУММАРНОЙ ОГИБАЮЩЕЙ
    // ============================================================
    
    double threshold = config_.envelope_threshold;
    if (threshold < 1e-6) threshold = 0.001;
    
    actual_duration_ = calculate_actual_duration_from_envelopes(threshold);
    
    // Не даём длительности стать слишком маленькой
    actual_duration_ = std::max(actual_duration_, 0.5);
    actual_duration_ = std::min(actual_duration_, global_duration);
    
    printf("   Actual duration (by envelope sum): %.2f s\n", actual_duration_);
    printf("   Threshold: %.6f\n", threshold);
    
    // ============================================================
    // 3. ОБРЕЗАЕМ ОГИБАЮЩИЕ ПО ФАКТИЧЕСКОЙ ДЛИТЕЛЬНОСТИ
    // ============================================================
    
    for (auto& pair : envelope_cache_) {
        EnvelopeData& data = pair.second;
        if (data.duration > actual_duration_) {
            data.duration = actual_duration_;
            size_t new_size = static_cast<size_t>(actual_duration_ / data.step) + 1;
            if (new_size < data.envelope.size()) {
                data.envelope.resize(new_size);
            }
        }
    }
    
    // ============================================================
    // 4. ГЕНЕРАЦИЯ СИГНАЛА
    // ============================================================
    
    size_t total_samples_os = static_cast<size_t>(actual_duration_ * 
                                                   config_.sample_rate * 
                                                   config_.oversample_factor);
    if (total_samples_os < 1) total_samples_os = 1;
    
    audio_oversampled_.assign(total_samples_os, 0.0);
    
    printf("\n   Synthesizing signal (sine modes, %zu modes)...\n", modes_.size());
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        auto it = envelope_cache_.find(m);
        if (it != envelope_cache_.end()) {
            add_mode(modes_[m], it->second, audio_oversampled_);
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
    // ПРИМЕНЕНИЕ ЭФФЕКТОВ
    // ============================================================
    
    printf("\n   Applying effects...\n");
    
    if (config_.reverb_mix > 0.001) {
        printf("     - Reverb (mix=%.2f)\n", config_.reverb_mix);
        add_reverb();
    }
    
    // Диагностика: после эффектов
    double max_effects = 0.0;
    size_t check_effects = std::min(size_t(1000), audio_.size());
    for (size_t i = 0; i < check_effects; ++i) {
        max_effects = std::max(max_effects, std::abs(audio_[i]));
    }
    printf("   audio_ after effects size = %zu, max amp = %.6f\n", 
           audio_.size(), max_effects);
    
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
        std::cerr << "Error: no audio to save\n";
        return false;
    }
    if (!Audio::writeMonoPcm16(filename, audio_, config_.sample_rate)) {
        std::cerr << "Error: failed to write valid PCM16 WAV file " << filename << "\n";
        return false;
    }
    const auto dataSize = audio_.size() * sizeof(std::int16_t);
    printf("   Saved: %s\n", filename.c_str());
    printf("   Size: %zu KB\n", dataSize / 1024u);
    return true;
}

// =========================================================================
// АНАЛИЗ СПЕКТРА
// =========================================================================

void BellSynthesizer::analyze_spectrum(const std::vector<double>& frequencies) {
    if (audio_.empty()) return;
    
    size_t n = audio_.size();
    
    printf("\n--- SPECTRUM ---\n");
    printf("   Target frequencies (Hz):\n");
    printf("   ------------------------\n");
    
    double time_norm = 1.0 / n;
    const double dt = 1.0 / config_.sample_rate;
    
    for (double freq : frequencies) {
        if (freq > config_.max_frequency) continue;
        
        double real = 0.0, imag = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double t = i * dt;
            double phase = 2.0 * bell::kPi * freq * t;
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

void BellSynthesizer::analyze_fft(int fft_size) {
    if (audio_.empty()) {
        std::cerr << "Error: no audio for FFT analysis\n";
        return;
    }

    const size_t max_fft_size = 65536;
    const size_t target = (fft_size <= 0)
        ? std::min(audio_.size(), max_fft_size)
        : std::min(static_cast<size_t>(fft_size), std::min(audio_.size(), max_fft_size));

    size_t n = 1;
    while (n < target) n <<= 1;
    if (n < 2) {
        std::cerr << "Error: FFT size too small (n < 2)\n";
        return;
    }

    printf("   Performing FFT (size: %zu, resolution: %.1f Hz)...\n",
           n, static_cast<double>(config_.sample_rate) / static_cast<double>(n));

    // In-place radix-2 Cooley-Tukey FFT. The project only needs a forward
    // real-signal spectrum, so a full external FFT package is unnecessary.
    struct Complex { double re, im; };
    std::vector<Complex> data(n, {0.0, 0.0});
    for (size_t i = 0; i < target; ++i) {
        const double window = 0.5 * (1.0 - std::cos(2.0 * bell::kPi * i / static_cast<double>(n - 1)));
        data[i].re = audio_[i] * window;
    }

    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        const double angle = -2.0 * bell::kPi / static_cast<double>(len);
        const double wlen_re = std::cos(angle);
        const double wlen_im = std::sin(angle);
        for (size_t i = 0; i < n; i += len) {
            double w_re = 1.0, w_im = 0.0;
            const size_t half = len >> 1;
            for (size_t j = 0; j < half; ++j) {
                const auto u = data[i + j];
                const auto v0 = data[i + j + half];
                const double v_re = v0.re * w_re - v0.im * w_im;
                const double v_im = v0.re * w_im + v0.im * w_re;
                data[i + j] = {u.re + v_re, u.im + v_im};
                data[i + j + half] = {u.re - v_re, u.im - v_im};
                const double next_re = w_re * wlen_re - w_im * wlen_im;
                w_im = w_re * wlen_im + w_im * wlen_re;
                w_re = next_re;
            }
        }
    }

    fft_frequencies_.clear();
    fft_magnitudes_.clear();
    const double window_compensation = 2.0;
    for (size_t i = 0; i < n / 2; ++i) {
        fft_frequencies_.push_back(static_cast<double>(i) * config_.sample_rate / static_cast<double>(n));
        const double mag = std::hypot(data[i].re, data[i].im) / static_cast<double>(n) * window_compensation;
        fft_magnitudes_.push_back(mag);
    }
}

void BellSynthesizer::print_fft_peaks(const std::vector<double>& target_frequencies) {
    if (fft_frequencies_.empty()) {
        std::cerr << "Please run analyze_fft() first\n";
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
        std::cerr << "Please run analyze_fft() first\n";
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed to create " << filename << "\n";
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
        std::cerr << "Error: envelopes not computed yet\n";
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed to create " << filename << "\n";
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