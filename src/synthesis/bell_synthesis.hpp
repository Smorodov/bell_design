// =============================================================================
// bell_synthesis.hpp - ПАРАМЕТРЫ СИНТЕЗА ЗВУКА КОЛОКОЛА
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <cmath>



namespace BellOptimizer {

struct Mode {
    double frequency;
    double amplitude;
    double decay_time;
    double phase;
    double attack_time;
    
    Mode() : frequency(0.0), amplitude(0.0), decay_time(0.0), 
             phase(0.0), attack_time(0.005) {}
};

struct GeometryParams {
    double height;
    double r_max;
    double r_min;
    double thickness;
    double crown_height = 0.0;
    double crown_radius = 0.0;
    double crown_thickness = 0.0;
    
    std::vector<double> profile_radius;
    std::vector<double> profile_thickness;
    bool has_profile = false;
    
    double get_r_max() const { return r_max; }
    double get_r_min() const { return r_min; }
    double get_radius_at(double height_pos) const;
    double get_thickness_at(double height_pos) const;
    bool is_monotonic() const;
    double get_max_slope() const;
};

struct SynthesisConfig {
    // Основные параметры
    double duration = 5.0;
    double fade_in = 0.003;
    double fade_out = 0.3;
    
    int sample_rate = 48000;
    int oversample_factor = 4;
    double max_frequency = 20000.0;
    
    double loss_factor = 0.0005;
    double tail_factor = 1.2;
    double envelope_table_step = 0.001;
    
    // Final output level. Modal ratios are fixed by the physical excitation
    // model; this parameter is only a global digital output scale.
    double target_peak = 0.95;
    
    // Параметры обрезки тишины
    double envelope_threshold = 0.00001;
    bool use_tail_trimming = true;
    
    // Реверберация
    double reverb_mix = 0.05;
    double reverb_feedback = 0.3;
    
    double target_fundamental = 440.0;
    std::string output_filename = "bell.wav";
    
    SynthesisConfig() = default;
};

class BellSynthesizer {
public:
    BellSynthesizer(const SynthesisConfig& config);
    ~BellSynthesizer() = default;
    
    // Основные методы
    void generate_modes(const GeometryParams& geo, const std::vector<Mode>& base_modes);
    void synthesize();
    
    // Сохранение
    bool save_wav(const std::string& filename);
    
    // Анализ
    void analyze_spectrum(const std::vector<double>& frequencies);
    void print_spectrum(const std::vector<Mode>& modes) const;
    void analyze_fft(int fft_size = 4096);
    void print_fft_peaks(const std::vector<double>& target_frequencies);
    void save_spectrum_csv(const std::string& filename, const std::vector<double>& target_frequencies);
    void save_all_envelopes_csv(const std::string& filename);
    
    // Геттеры
    const std::vector<double>& get_audio() const { return audio_; }
    const std::vector<Mode>& get_modes() const { return modes_; }
    int get_sample_rate() const { return config_.sample_rate; }
    double get_actual_duration() const { return actual_duration_; }
    
private:
    struct EnvelopeData {
        std::vector<double> envelope;
        double step;
        double duration;
    };
    
    // Внутренние методы
    double calculate_natural_duration() const;
    void precompute_envelopes(double duration);
    void downsample_audio();
    
    void add_mode(const Mode& mode, const EnvelopeData& envelope_data,
                  std::vector<double>& buffer);
    
    double calculate_decay_envelope(double time, const Mode& mode) const;
    double get_envelope_value(const EnvelopeData& data, double time) const;
    
    // ✅ НОВЫЙ МЕТОД ДЛЯ РАСЧЕТА ДЛИТЕЛЬНОСТИ ПО СУММАРНОЙ ОГИБАЮЩЕЙ
    double calculate_actual_duration_from_envelopes(double threshold);
    
    // Эффекты
    void add_reverb();
    void apply_fade(double fade_in_time, double fade_out_time);
    void normalize(double peak_level);
    void trim_silence(double threshold);
    
    void write_wav_header(std::ofstream& file, int data_size) const;
    
    SynthesisConfig config_;
    std::vector<Mode> modes_;
    std::vector<double> audio_;
    std::vector<double> audio_oversampled_;
    GeometryParams geometry_;
    std::map<int, EnvelopeData> envelope_cache_;
    
    double actual_duration_ = 0.0;
    double peak_amplitude_ = 0.0;
    double rms_amplitude_ = 0.0;
    
    std::vector<double> fft_frequencies_;
    std::vector<double> fft_magnitudes_;
};

} // namespace BellOptimizer