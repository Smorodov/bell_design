// =============================================================================
// bell_synthesis.cpp - РЕАЛИЗАЦИЯ СИНТЕЗА ЗВУКА КОЛОКОЛА
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

double GeometryParams::get_radius_at(double height_pos) const {
    if (!has_profile || profile_radius.empty()) {
        double t = std::max(0.0, std::min(1.0, height_pos));
        return r_min + (r_max - r_min) * t;
    }
    
    double pos = height_pos * (profile_radius.size() - 1);
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

BellSynthesizer::BellSynthesizer(const SynthesisConfig& config)
    : config_(config) {
    if (config_.duration > 0.0) {
        audio_.reserve(static_cast<size_t>(config_.duration * config_.sample_rate));
    }
}

void BellSynthesizer::generate_modes(const GeometryParams& geo, 
                                     const std::vector<Mode>& base_modes) {
    geometry_ = geo;
    modes_.clear();
    
    if (!base_modes.empty()) {
        modes_ = base_modes;
    } else {
        fmt::print(stderr, "Ошибка: нет мод для синтеза\n");
        return;
    }
    
    fmt::print("   Мод для синтеза: {}\n", modes_.size());
}

void BellSynthesizer::generate_modes_from_params(const GeometryParams& geo) {
    fmt::print(stderr, "Warning: generate_modes_from_params() вызван, но моды должны браться только из .dat\n");
}

double BellSynthesizer::calculate_decay_envelope(double time, const Mode& mode) const {
    double t60 = mode.decay_time;
    if (t60 < 0.001) t60 = 0.001;
    
    double attack_time = mode.attack_time;
    if (attack_time < 0.001) attack_time = 0.001;
    
    double attack = 1.0 - std::exp(-time / attack_time);
    
    double amplitude = mode.amplitude * attack;
    double nonlinear_factor = 1.0 + config_.nonlinear_damping * amplitude * amplitude;
    double effective_t60 = t60 / nonlinear_factor;
    
    double decay = std::exp(-3.0 * time / effective_t60);
    
    double freq_factor = mode.frequency / 1000.0;
    double correction = 1.0 + 0.1 * std::exp(-freq_factor);
    
    double t_rel = time / effective_t60;
    
    if (t_rel < 0.15) {
        double freq_boost = 1.0 + config_.boost_freq_factor * freq_factor;
        double boost_factor = 1.0 + config_.boost_amplitude * freq_boost * (1.0 - t_rel / 0.15);
        decay = decay * boost_factor;
    }
    
    if (t_rel > 0.7) {
        double slowdown = 1.0 + config_.slowdown_factor * ((t_rel - 0.7) / 0.3);
        decay = std::pow(decay, 1.0 / (slowdown * correction));
    }
    
    double envelope = decay * attack;
    
    if (envelope < 1e-12) envelope = 0.0;
    
    return envelope;
}

double BellSynthesizer::get_envelope_value(const EnvelopeData& data, double time) const {
    if (data.envelope.empty()) return 0.0;
    if (time < 0.0) return data.envelope[0];
    if (time >= data.duration) return 0.0;
    
    double idx_double = time / data.step;
    size_t idx = static_cast<size_t>(idx_double);
    double frac = idx_double - idx;
    
    if (idx >= data.envelope.size() - 1) return data.envelope.back();
    
    return data.envelope[idx] * (1.0 - frac) + data.envelope[idx + 1] * frac;
}

std::vector<double> BellSynthesizer::scale_envelope_to_duration(
    const std::vector<double>& original_envelope,
    double original_duration,
    double target_duration,
    const Mode& mode) {
    
    if (original_envelope.empty()) return original_envelope;
    if (target_duration <= 0.0) return original_envelope;
    
    if (std::abs(original_duration - target_duration) < 1e-6) {
        return original_envelope;
    }
    
    // Находим естественную длительность (где огибающая падает ниже порога)
    double peak = mode.amplitude;
    double noise_floor = 1e-8;
    double relative_threshold = config_.envelope_threshold;
    double threshold = std::max(noise_floor, peak * relative_threshold);
    if (threshold < 1e-6) threshold = 1e-6;
    
    bool was_above_threshold = false;
    size_t natural_end = original_envelope.size();
    
    for (size_t i = 0; i < original_envelope.size(); ++i) {
        if (original_envelope[i] > threshold) {
            was_above_threshold = true;
        }
        if (was_above_threshold && original_envelope[i] < threshold) {
            natural_end = i;
            break;
        }
    }
    
    if (natural_end == original_envelope.size()) {
        natural_end = (size_t)(original_envelope.size() * 0.99);
    }
    if (natural_end < original_envelope.size() * 0.1) {
        natural_end = (size_t)(original_envelope.size() * 0.1);
    }
    
    // Масштабируем только активную часть (0 ... natural_end)
    double original_active_duration = (double)natural_end / (original_envelope.size() - 1) * original_duration;
    double scale = target_duration / original_active_duration;
    
    if (original_active_duration < 1e-12) {
        return original_envelope;
    }
    
    std::vector<double> scaled_envelope;
    scaled_envelope.resize(original_envelope.size());
    
    for (size_t i = 0; i < original_envelope.size(); ++i) {
        double t_norm = (double)i / (original_envelope.size() - 1);
        double t_scaled = t_norm / scale;
        
        if (t_scaled >= 1.0) {
            scaled_envelope[i] = 0.0;  // Хвост обнуляем
        } else {
            double pos = t_scaled * (original_envelope.size() - 1);
            size_t idx = (size_t)pos;
            double frac = pos - idx;
            
            if (idx >= original_envelope.size() - 1) {
                scaled_envelope[i] = 0.0;
            } else {
                scaled_envelope[i] = original_envelope[idx] * (1.0 - frac) + 
                                     original_envelope[idx + 1] * frac;
            }
        }
    }
    
    // Энергетическая компенсация НЕ ПРИМЕНЯЕТСЯ
    // Нормализация на выходе выровняет громкость всех файлов
    
    return scaled_envelope;
}

void BellSynthesizer::precompute_envelopes() {
    envelope_cache_.clear();
    
    bool use_natural_duration = (config_.duration < 0.0);
    
    double max_t60 = 0.0;
    for (const auto& mode : modes_) {
        if (mode.decay_time > max_t60) max_t60 = mode.decay_time;
    }
    
    double natural_duration = max_t60 * config_.tail_factor;
    if (natural_duration > 30.0) natural_duration = 30.0;
    if (natural_duration < 1.0) natural_duration = 3.0;
    
    if (use_natural_duration || config_.duration <= 0.0) {
        actual_duration_ = natural_duration;
    } else {
        actual_duration_ = config_.duration;
        if (actual_duration_ < 0.1) actual_duration_ = 0.1;
        if (actual_duration_ > 30.0) actual_duration_ = 30.0;
    }
    
    double table_step = config_.envelope_table_step;
    if (table_step < 0.0001) table_step = 0.0001;
    
    bool need_scaling = !use_natural_duration && 
                        std::abs(natural_duration - actual_duration_) > 1e-6;
    
    fmt::print("   Расчёт огибающих (шаг {:.1f} мс)...\n", table_step * 1000.0);
    if (need_scaling) {
        fmt::print("   Естественная длительность: {:.2f} с\n", natural_duration);
        fmt::print("   Целевая длительность: {:.2f} с (масштабирование)\n", actual_duration_);
    }
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        EnvelopeData data;
        data.step = table_step;
        data.duration = actual_duration_;
        
        if (need_scaling) {
            size_t natural_size = static_cast<size_t>(natural_duration / table_step) + 1;
            if (natural_size < 10) natural_size = 10;
            
            std::vector<double> natural_envelope(natural_size);
            for (size_t i = 0; i < natural_size; ++i) {
                double t = i * table_step;
                natural_envelope[i] = calculate_decay_envelope(t, modes_[m]);
            }
            
            data.envelope = scale_envelope_to_duration(
                natural_envelope, natural_duration, actual_duration_, modes_[m]);
        } else {
            size_t table_size = static_cast<size_t>(actual_duration_ / table_step) + 1;
            if (table_size < 10) table_size = 10;
            data.envelope.resize(table_size);
            for (size_t i = 0; i < table_size; ++i) {
                double t = i * table_step;
                data.envelope[i] = calculate_decay_envelope(t, modes_[m]);
            }
        }
        
        envelope_cache_[m] = std::move(data);
    }
    
    fmt::print("   Огибающие рассчитаны: {} мод\n", modes_.size());
    fmt::print("   Длительность: {:.2f} с", actual_duration_);
    if (use_natural_duration) {
        fmt::print(" (естественная)");
    } else if (need_scaling) {
        fmt::print(" (масштабирована из {:.2f} с)", natural_duration);
    } else {
        fmt::print(" (фиксированная)");
    }
    fmt::print("\n");
}

void BellSynthesizer::add_mode_to_audio(const Mode& mode, 
                                        const EnvelopeData& envelope_data,
                                        std::vector<double>& buffer) {
    if (envelope_data.envelope.empty()) return;
    
    double dt = 1.0 / (config_.sample_rate * config_.oversample_factor);
    size_t max_samples = std::min(buffer.size(), 
                                  static_cast<size_t>(envelope_data.duration / dt));
    double amp = mode.amplitude;
    double freq = mode.frequency;
    
    double nyquist = 0.5 / dt;
    if (freq >= nyquist) {
        fmt::print(stderr, "Warning: частота {:.1f} Гц превышает Найквиста ({:.1f} Гц), пропускаем моду\n", 
                   freq, nyquist);
        return;
    }
    
    for (size_t i = 0; i < max_samples; ++i) {
        double t = i * dt;
        double envelope = get_envelope_value(envelope_data, t);
        double oscillation = std::sin(2.0 * M_PI * mode.frequency * t + mode.phase);
        buffer[i] += amp * envelope * oscillation;
    }
}

void BellSynthesizer::synthesize() {
    if (modes_.empty()) {
        fmt::print(stderr, "Ошибка: нет мод для синтеза\n");
        return;
    }
    
    precompute_envelopes();
    
    size_t total_samples_os = static_cast<size_t>(actual_duration_ * 
                                                   config_.sample_rate * 
                                                   config_.oversample_factor);
    if (total_samples_os < 1) total_samples_os = 1;
    
    audio_oversampled_.assign(total_samples_os, 0.0);
    
    fmt::print("   Синтез звука ({:.2f} с, {} Гц, x{} oversample)...\n", 
               actual_duration_, config_.sample_rate, config_.oversample_factor);
    
    for (size_t m = 0; m < modes_.size(); ++m) {
        auto it = envelope_cache_.find(m);
        if (it != envelope_cache_.end()) {
            add_mode_to_audio(modes_[m], it->second, audio_oversampled_);
        }
    }
    
    size_t total_samples = static_cast<size_t>(actual_duration_ * config_.sample_rate);
    if (total_samples < 1) total_samples = 1;
    
    audio_.assign(total_samples, 0.0);
    
    for (size_t i = 0; i < total_samples; ++i) {
        double sum = 0.0;
        size_t start = i * config_.oversample_factor;
        size_t end = std::min(start + config_.oversample_factor, total_samples_os);
        
        for (size_t j = start; j < end; ++j) {
            sum += audio_oversampled_[j];
        }
        audio_[i] = sum / (end - start);
    }
    
    if (config_.use_acoustic_model) {
        apply_acoustic_coupling();
        if (config_.use_helmholtz_resonance) {
            apply_helmholtz_resonance();
        }
        apply_radiation_filter();
    }
    
    if (config_.use_mount_model) {
        apply_mount_filter();
    }
    
    fmt::print("   Синтез завершён\n");
}

void BellSynthesizer::apply_acoustic_coupling() {
    if (!config_.use_acoustic_coupling) return;
    
    size_t total_samples = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    for (const auto& mode : modes_) {
        if (mode.frequency > config_.max_frequency) continue;
        
        double amp = mode.amplitude * 0.05;
        double freq = mode.frequency;
        double damping = 0.01;
        
        for (size_t i = 0; i < total_samples; ++i) {
            double t = i * dt;
            double envelope = std::exp(-damping * t);
            double osc = std::sin(2.0 * M_PI * freq * t + mode.phase);
            audio_[i] += amp * envelope * osc;
        }
    }
}

void BellSynthesizer::apply_helmholtz_resonance() {
    if (!config_.use_helmholtz_resonance) return;
    
    size_t total_samples = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    double V = M_PI * geometry_.r_max * geometry_.r_max * geometry_.height;
    double A = M_PI * geometry_.r_max * geometry_.r_max;
    double L = geometry_.height * 0.1;
    double c = 343.0;
    
    double f_h = (c / (2.0 * M_PI)) * std::sqrt(A / (V * L));
    f_h = std::max(f_h, 50.0);
    f_h = std::min(f_h, 500.0);
    
    double amp = config_.helmholtz_amplitude;
    double damping = config_.helmholtz_damping;
    
    for (size_t i = 0; i < total_samples; ++i) {
        double t = i * dt;
        double envelope = std::exp(-damping * t);
        double osc = std::sin(2.0 * M_PI * f_h * t);
        audio_[i] += amp * envelope * osc;
    }
}

void BellSynthesizer::apply_radiation_filter() {
    if (!config_.use_acoustic_model) return;
    
    double efficiency = config_.radiation_efficiency;
    double damping = config_.radiation_damping_factor;
    
    size_t total_samples = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    for (size_t i = 0; i < total_samples; ++i) {
        double t = i * dt;
        double filter = 1.0 + efficiency * std::exp(-damping * t);
        audio_[i] *= filter;
    }
}

void BellSynthesizer::apply_mount_filter() {
    if (!config_.use_mount_model) return;
    
    size_t total_samples = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    double damping = config_.mount_damping;
    
    double fc = config_.target_fundamental * 0.5;
    double alpha = 1.0 / (1.0 + 2.0 * M_PI * fc * dt);
    
    double prev = 0.0;
    for (size_t i = 0; i < total_samples; ++i) {
        double sample = audio_[i];
        double filtered = alpha * sample + (1.0 - alpha) * prev;
        audio_[i] = filtered * (1.0 - damping * (1.0 - filtered / (sample + 1e-10)));
        prev = filtered;
    }
}

void BellSynthesizer::add_reverb() {
    if (config_.reverb_mix < 0.001) return;
    
    size_t total_samples = audio_.size();
    
    int delay_samples = static_cast<int>(config_.reverb_delay * config_.sample_rate);
    double feedback = config_.reverb_feedback;
    double mix = config_.reverb_mix;
    
    if (delay_samples >= static_cast<int>(total_samples)) {
        delay_samples = static_cast<int>(total_samples * 0.05);
    }
    
    std::vector<double> reverb(total_samples, 0.0);
    
    for (size_t i = delay_samples; i < total_samples; ++i) {
        reverb[i] = audio_[i - delay_samples] * 0.7 + 
                    reverb[i - delay_samples] * feedback * 0.3;
    }
    
    for (size_t i = 0; i < total_samples; ++i) {
        audio_[i] = audio_[i] * (1.0 - mix) + reverb[i] * mix;
    }
}

void BellSynthesizer::apply_fade(double fade_in_time, double fade_out_time) {
    size_t total_samples = audio_.size();
    
    size_t fade_in_samples = static_cast<size_t>(fade_in_time * config_.sample_rate);
    size_t fade_out_samples = static_cast<size_t>(fade_out_time * config_.sample_rate);
    
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

bool BellSynthesizer::save_wav(const std::string& filename) {
    if (audio_.empty()) {
        fmt::print(stderr, "Ошибка: нет аудио для сохранения\n");
        return false;
    }
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        fmt::print(stderr, "Ошибка: не удалось создать файл {}\n", filename);
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
        fmt::print(stderr, "Ошибка: запись в файл не удалась\n");
        file.close();
        return false;
    }
    
    file.close();
    
    fmt::print("   Сохранён: {}\n", filename);
    fmt::print("   Размер: {} KB\n", data_size / 1024);
    
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

void BellSynthesizer::analyze_spectrum(const std::vector<double>& frequencies) {
    if (audio_.empty()) return;
    
    size_t n = audio_.size();
    double dt = 1.0 / config_.sample_rate;
    
    fmt::print("\n🎵 СПЕКТР СИНТЕЗИРОВАННОГО ЗВУКА:\n");
    fmt::print("   {:-<60}\n", "");
    fmt::print("   Основные частоты (Гц):\n");
    fmt::print("   {:-<60}\n", "");
    
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
            fmt::print("   {:.1f} Гц : {:.3f} мВ\n", freq, amp * 1000.0);
        }
    }
    
    fmt::print("   {:-<60}\n", "");
}

void BellSynthesizer::print_spectrum(const std::vector<Mode>& modes) const {
    if (modes.empty()) return;
    
    fmt::print("\n🎵 СПЕКТР МОД:\n");
    fmt::print("   {:-<60}\n", "");
    fmt::print("   № | Частота(Гц) | T60(с) | Атака(с) | Амплитуда\n");
    fmt::print("   {:-<60}\n", "");
    
    for (size_t i = 0; i < modes.size() && i < 25; ++i) {
        const auto& m = modes[i];
        fmt::print("   {:2d} | {:10.1f} | {:6.1f} | {:7.3f} | {:8.6f}\n",
                   i+1, m.frequency, m.decay_time, m.attack_time, m.amplitude);
    }
    fmt::print("   {:-<60}\n", "");
}

void BellSynthesizer::analyze_fft(int fft_size) {
    if (audio_.empty()) {
        fmt::print(stderr, "Ошибка: нет аудио для FFT анализа\n");
        return;
    }

#ifndef USE_FFTW
    fmt::print(stderr, "Ошибка: FFTW не доступен. Пересоберите проект с USE_FFTW\n");
    return;
#endif

    size_t max_fft_size = 65536;
    size_t target = (fft_size <= 0) ? std::min(audio_.size(), max_fft_size) 
                                    : std::min(static_cast<size_t>(fft_size), 
                                               std::min(audio_.size(), max_fft_size));
    
    size_t n = 1;
    while (n < target) n <<= 1;

    if (n < 2) {
        fmt::print(stderr, "Ошибка: размер FFT слишком мал (n < 2)\n");
        return;
    }

    fmt::print("   Выполняется FFT (размер: {}, разрешение: {:.1f} Гц)...\n", 
               n, (config_.sample_rate / n));

#ifdef USE_FFTW
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    
    if (!in || !out) {
        fmt::print(stderr, "Ошибка: не удалось выделить память для FFT\n");
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
        double freq = i / (n * dt);
        double mag = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) / n * window_compensation;
        fft_frequencies_.push_back(freq);
        fft_magnitudes_.push_back(mag);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
#endif
}

void BellSynthesizer::print_fft_peaks(const std::vector<double>& target_frequencies) {
    if (fft_frequencies_.empty()) {
        fmt::print(stderr, "Сначала выполните analyze_fft()\n");
        return;
    }

    fmt::print("\n📊 СПЕКТР (ТОЛЬКО ЗАДАННЫЕ ЧАСТОТЫ):\n");
    fmt::print("   {:-<60}\n", "");
    fmt::print("   Частота(Гц) | Амплитуда\n");
    fmt::print("   {:-<60}\n", "");

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
            fmt::print("   {:.2f} Гц -> {:.6f} (найдено {:.2f} Гц)\n", 
                       target_freq, best_mag, best_freq);
        } else {
            fmt::print("   {:.2f} Гц -> НЕ НАЙДЕНА\n", target_freq);
        }
    }
    
    fmt::print("   {:-<60}\n", "");
}

void BellSynthesizer::save_spectrum_csv(const std::string& filename, 
                                        const std::vector<double>& target_frequencies) {
    if (fft_frequencies_.empty()) {
        fmt::print(stderr, "Сначала выполните analyze_fft()\n");
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        fmt::print(stderr, "Ошибка: не удалось создать {}\n", filename);
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
    
    fmt::print("   Спектр сохранён: {} ({} из {} частот найдено)\n", 
               filename, found_count, target_frequencies.size());
}

void BellSynthesizer::save_all_envelopes_csv(const std::string& filename) {
    if (envelope_cache_.empty()) {
        fmt::print(stderr, "Ошибка: огибающие ещё не рассчитаны\n");
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        fmt::print(stderr, "Ошибка: не удалось создать {}\n", filename);
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
    fmt::print("   Все огибающие сохранены: {}\n", filename);
}

}