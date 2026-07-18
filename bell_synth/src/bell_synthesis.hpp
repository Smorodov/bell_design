// =============================================================================
// bell_synthesis.hpp - ПАРАМЕТРЫ СИНТЕЗА ЗВУКА
// =============================================================================
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <map>
#include <complex>

#include <Eigen/Dense>

#include "bell_constants.hpp"

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
    double duration = 5.0;
    double fade_in = BellConstants::ATTACK_DURATION;
    double fade_out = 0.2;
    
    int sample_rate = 48000;
    int oversample_factor = 4;
    
    int max_overtone_count = 32;
    double min_overtone_amplitude = BellConstants::MIN_AMPLITUDE;
    double max_frequency = BellConstants::MAX_FREQUENCY;
    
    double t60_base = 8.0;
    double alpha_damping = 0.0001;
    double beta_damping = 0.05;
    
    bool use_mount_model = true;
    double mount_damping = 0.2;
    double mount_freq_shift = 0.98;
    double mount_high_freq_attenuation = 0.8;
    
    bool use_acoustic_model = true;
    bool use_acoustic_coupling = true;
    bool use_helmholtz_resonance = true;
    double helmholtz_amplitude = 0.04;
    double helmholtz_damping = 0.3;
    double radiation_efficiency = 0.4;
    double radiation_damping_factor = 0.005;
    
    double reverb_mix = 0.12;
    double reverb_delay = 0.1;
    double reverb_feedback = 0.5;
    
    double target_peak = BellConstants::CLIP_LEVEL;
    double target_rms = 0.35;
    bool use_rms_normalization = true;
    bool soft_clip = true;
    
    double impact_position = 0.72;
    double impact_duration = 0.0003;
    double impact_force = 1.0;
    
    double target_fundamental = 440.0;
    std::vector<double> target_intervals = {
        0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00
    };
    std::vector<double> target_amplitudes = {
        0.80, 0.40, 1.00, 0.60, 0.40, 0.30, 0.20, 0.15, 0.10, 0.08
    };
    
    double min_frequency = BellConstants::MIN_FREQUENCY;
    
    std::string output_filename = "orthodox_bell.wav";
    
    double boost_amplitude = 0.2;
    double boost_freq_factor = 0.5;
    double slowdown_factor = 0.15;
    
    double loss_factor = 2e-4;
    double nonlinear_damping = 0.5;
    
    double envelope_table_step = 0.001;
    
    double envelope_threshold = 1e-3;
    double tail_factor = 3.0;
    
    SynthesisConfig() = default;
};

class BellSynthesizer {
public:
    BellSynthesizer(const SynthesisConfig& config);
    ~BellSynthesizer() = default;
    
    void generate_modes(const GeometryParams& geo, const std::vector<Mode>& base_modes);
    void generate_modes_from_params(const GeometryParams& geo);
    
    void synthesize();
    void add_reverb();
    void apply_fade(double fade_in_time, double fade_out_time);
    void normalize(double peak_level);
    
    bool save_wav(const std::string& filename);
    void analyze_spectrum(const std::vector<double>& frequencies);
    void print_spectrum(const std::vector<Mode>& modes) const;
    
    void analyze_fft(int fft_size = 4096);
    void print_fft_peaks(const std::vector<double>& target_frequencies);
    void save_spectrum_csv(const std::string& filename, const std::vector<double>& target_frequencies);
    
    void save_all_envelopes_csv(const std::string& filename);
    
    const std::vector<double>& get_audio() const { return audio_; }
    const std::vector<Mode>& get_modes() const { return modes_; }
    int get_modes_count() const { return static_cast<int>(modes_.size()); }
    double get_optimal_impact_position() const { return config_.impact_position; }
    int get_sample_rate() const { return config_.sample_rate; }
    double get_actual_duration() const { return actual_duration_; }
    
private:
    struct EnvelopeData {
        std::vector<double> envelope;
        double step;
        double duration;
    };
    
    double calculate_decay_envelope(double time, const Mode& mode) const;
    double get_envelope_value(const EnvelopeData& data, double time) const;
    void add_mode_to_audio(const Mode& mode, const EnvelopeData& envelope_data, std::vector<double>& buffer);
    void precompute_envelopes();
    std::vector<double> scale_envelope_to_duration(
        const std::vector<double>& original_envelope,
        double original_duration,
        double target_duration,
        const Mode& mode);
    
    void apply_mount_filter();
    void apply_acoustic_coupling();
    void apply_helmholtz_resonance();
    void apply_radiation_filter();
    
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

}