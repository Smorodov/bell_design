/* ------------------------------------------------------------------
* Copyright (C) 2020 ewan xu<ewan_xu@outlook.com>
* Full implementation with Eigen - All features preserved
* All critical bugs fixed
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*      http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
* express or implied.
* See the License for the specific language governing permissions
* and limitations under the License.
* -------------------------------------------------------------------
*/

#ifndef LIBROSA_H_
#define LIBROSA_H_

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <vector>
#include <complex>
#include <iostream>
#include <random>
#include <cmath>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <limits>

namespace librosa {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr float EPSILON = 1e-10f;
constexpr float WINDOW_SUM_EPSILON = 1e-10f;
constexpr float DB_LIMIT = 80.0f;

namespace internal {
typedef Eigen::Matrix<float, 1, Eigen::Dynamic, Eigen::RowMajor> Vectorf;
typedef Eigen::Matrix<std::complex<float>, 1, Eigen::Dynamic, Eigen::RowMajor> Vectorcf;
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrixf;
typedef Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrixcf;
}

struct SpectrogramData {
    std::vector<std::vector<float>> magnitude;
    std::vector<std::vector<std::complex<float>>> phase;
    int n_frames;
    int n_freq;
    
    SpectrogramData() : n_frames(0), n_freq(0) {}
    SpectrogramData(int frames, int freq) : n_frames(frames), n_freq(freq) {
        magnitude.resize(frames, std::vector<float>(freq, 0.0f));
        phase.resize(frames, std::vector<std::complex<float>>(freq, 0.0f));
    }
    
    float& mag(int frame, int freq) { return magnitude[frame][freq]; }
    const float& mag(int frame, int freq) const { return magnitude[frame][freq]; }
    
    std::complex<float>& ph(int frame, int freq) { return phase[frame][freq]; }
    const std::complex<float>& ph(int frame, int freq) const { return phase[frame][freq]; }
    
    internal::Matrixf to_eigen_mag() const {
        internal::Matrixf result(n_frames, n_freq);
        for (int i = 0; i < n_frames; ++i) {
            for (int j = 0; j < n_freq; ++j) {
                result(i, j) = magnitude[i][j];
            }
        }
        return result;
    }
    
    internal::Matrixcf to_eigen_complex() const {
        internal::Matrixcf result(n_frames, n_freq);
        for (int i = 0; i < n_frames; ++i) {
            for (int j = 0; j < n_freq; ++j) {
                float mag = magnitude[i][j];
                float phase_angle = std::arg(phase[i][j]);
                result(i, j) = std::polar(mag, phase_angle);
            }
        }
        return result;
    }
};

namespace internal {

template<typename MatrixType>
static std::unique_ptr<MatrixType> create_matrix(int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        return std::make_unique<MatrixType>(0, 0);
    }
    
    if (rows > 1000000 || cols > 1000000) {
        throw std::runtime_error("Matrix dimensions too large");
    }
    
    try {
        auto matrix = std::make_unique<MatrixType>(rows, cols);
        matrix->setZero();
        return matrix;
    } catch (const std::bad_alloc& e) {
        throw std::runtime_error("Failed to allocate memory for matrix");
    }
}

static Vectorf to_eigen_vector(const std::vector<float>& vec) {
    Vectorf result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result(static_cast<Eigen::Index>(i)) = vec[i];
    }
    return result;
}

static std::vector<float> from_eigen_vector(const Vectorf& vec) {
    std::vector<float> result(vec.size());
    for (Eigen::Index i = 0; i < vec.size(); ++i) {
        result[static_cast<size_t>(i)] = vec(i);
    }
    return result;
}

static Vectorf hann_window(int n_fft) {
    if (n_fft <= 1) {
        throw std::invalid_argument("n_fft must be greater than 1");
    }
    
    Vectorf window(n_fft);
    for (int i = 0; i < n_fft; ++i) {
        window(i) = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (n_fft - 1)));
    }
    return window;
}

static Vectorf hamming_window(int n_fft) {
    if (n_fft <= 1) {
        throw std::invalid_argument("n_fft must be greater than 1");
    }
    
    Vectorf window(n_fft);
    for (int i = 0; i < n_fft; ++i) {
        window(i) = 0.54f - 0.46f * std::cos(2.0f * M_PI * i / (n_fft - 1));
    }
    return window;
}

static Vectorf blackman_window(int n_fft) {
    if (n_fft <= 1) {
        throw std::invalid_argument("n_fft must be greater than 1");
    }
    
    Vectorf window(n_fft);
    float a0 = 0.42f;
    float a1 = 0.5f;
    float a2 = 0.08f;
    for (int i = 0; i < n_fft; ++i) {
        window(i) = a0 - a1 * std::cos(2.0f * M_PI * i / (n_fft - 1)) + 
                    a2 * std::cos(4.0f * M_PI * i / (n_fft - 1));
    }
    return window;
}

static Vectorf generate_window(const std::string& win, int n_fft) {
    if (win == "hann" || win == "hanning") {
        return hann_window(n_fft);
    } else if (win == "hamming") {
        return hamming_window(n_fft);
    } else if (win == "blackman") {
        return blackman_window(n_fft);
    } else {
        throw std::invalid_argument("Unsupported window type: " + win);
    }
}

static Vectorf pad_reflect(const Vectorf& x, int left, int right) {
    Vectorf x_paded = Vectorf::Constant(left + x.size() + right, 0.0f);
    x_paded.segment(left, x.size()) = x;
    
    if (x.size() == 0) {
        return x_paded;
    }
    
    for (int i = 0; i < left; ++i) {
        int src_idx = left - i;
        if (src_idx >= 0 && src_idx < static_cast<int>(x.size())) {
            x_paded[i] = x[src_idx];
        } else {
            x_paded[i] = x[0];
        }
    }
    
    for (int i = 0; i < right; ++i) {
        int src_idx = static_cast<int>(x.size()) - 2 - i;
        if (src_idx >= 0 && src_idx < static_cast<int>(x.size())) {
            x_paded[left + x.size() + i] = x[src_idx];
        } else {
            x_paded[left + x.size() + i] = x[x.size() - 1];
        }
    }
    
    return x_paded;
}

static Vectorf pad_symmetric(const Vectorf& x, int left, int right) {
    Vectorf x_paded = Vectorf::Constant(left + x.size() + right, 0.0f);
    x_paded.segment(left, x.size()) = x;
    
    if (x.size() == 0) {
        return x_paded;
    }
    
    for (int i = 0; i < left; ++i) {
        int src_idx = left - i - 1;
        if (src_idx >= 0 && src_idx < static_cast<int>(x.size())) {
            x_paded[i] = x[src_idx];
        } else {
            x_paded[i] = x[0];
        }
    }
    
    for (int i = 0; i < right; ++i) {
        int src_idx = static_cast<int>(x.size()) - 1 - i;
        if (src_idx >= 0 && src_idx < static_cast<int>(x.size())) {
            x_paded[left + x.size() + i] = x[src_idx];
        } else {
            x_paded[left + x.size() + i] = x[x.size() - 1];
        }
    }
    
    return x_paded;
}

static Vectorf pad_edge(const Vectorf& x, int left, int right) {
    Vectorf x_paded = Vectorf::Constant(left + x.size() + right, 0.0f);
    x_paded.segment(left, x.size()) = x;
    
    if (x.size() == 0) {
        return x_paded;
    }
    
    for (int i = 0; i < left; ++i) {
        x_paded[i] = x[0];
    }
    
    for (int i = 0; i < right; ++i) {
        x_paded[left + x.size() + i] = x[x.size() - 1];
    }
    
    return x_paded;
}

static Vectorf pad(const Vectorf& x, int left, int right, const std::string& mode, float value) {
    if (left < 0 || right < 0) {
        throw std::invalid_argument("Padding sizes must be non-negative");
    }
    
    if (left > 1000000 || right > 1000000) {
        throw std::invalid_argument("Padding sizes too large");
    }
    
    if (mode == "reflect") {
        return pad_reflect(x, left, right);
    } else if (mode == "symmetric") {
        return pad_symmetric(x, left, right);
    } else if (mode == "edge") {
        return pad_edge(x, left, right);
    } else if (mode == "constant") {
        Vectorf x_paded = Vectorf::Constant(left + x.size() + right, value);
        x_paded.segment(left, x.size()) = x;
        return x_paded;
    } else {
        throw std::invalid_argument("Unsupported padding mode: " + mode);
    }
}

static void validate_mel_params(int sr, int n_fft, int n_mels, int fmin, int fmax) {
    if (sr <= 0) {
        throw std::invalid_argument("Sample rate must be positive");
    }
    if (n_fft <= 0 || n_fft % 2 != 0) {
        throw std::invalid_argument("n_fft must be positive even number");
    }
    if (n_mels <= 0) {
        throw std::invalid_argument("Number of mel bands must be positive");
    }
    if (fmin < 0) {
        throw std::invalid_argument("fmin must be non-negative");
    }
    if (fmax < 0 || fmax > sr / 2) {
        throw std::invalid_argument("fmax must be between 0 and sr/2");
    }
    if (fmin >= fmax) {
        throw std::invalid_argument("fmin must be less than fmax");
    }
}

static Matrixcf stft_impl(const Vectorf& x, int n_fft, int n_hop,
                         const Vectorf& window, bool center,
                         const std::string& mode) {
    if (n_fft <= 0 || n_hop <= 0) {
        throw std::invalid_argument("n_fft and n_hop must be positive");
    }
    
    if (window.size() != n_fft) {
        throw std::invalid_argument("Window size must match n_fft");
    }
    
    if (x.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("Input signal too large");
    }
    
    int pad_len = center ? n_fft / 2 : 0;
    Vectorf x_paded = pad(x, pad_len, pad_len, mode, 0.0f);
    
    int n_f = n_fft / 2 + 1;
    
    size_t signal_size = center ? x_paded.size() : x.size();
    int n_frames = 0;
    
    if (center) {
        if (signal_size >= static_cast<size_t>(n_fft)) {
            n_frames = static_cast<int>((signal_size - n_fft) / n_hop) + 1;
        } else {
            n_frames = 0;
        }
    } else {
        if (x.size() >= static_cast<size_t>(n_fft)) {
            n_frames = static_cast<int>((x.size() - n_fft) / n_hop) + 1;
        } else {
            n_frames = 0;
        }
    }
    
    if (n_frames == 0) {
        return Matrixcf(0, n_f);
    }
    
    if (n_frames > 1000000) {
        throw std::invalid_argument("Too many frames");
    }
    
    auto X = create_matrix<Matrixcf>(n_frames, n_f);
    Eigen::FFT<float> fft;
    
    Vectorf x_frame(n_fft);
    Vectorcf X_full(n_fft);
    
    for (int i = 0; i < n_frames; ++i) {
        int start = i * n_hop;
        x_frame.setZero();
        
        int copy_len = std::min(n_fft, static_cast<int>(x_paded.size() - start));
        if (copy_len > 0 && start + copy_len <= static_cast<int>(x_paded.size())) {
            x_frame.head(copy_len) = x_paded.segment(start, copy_len);
        }
        
        x_frame = window.array() * x_frame.array();
        X_full = fft.fwd(x_frame);
        (*X).row(i) = X_full.head(n_f);
    }
    
    return std::move(*X);
}

static Vectorf istft_with_phases(const Matrixcf& X, int n_fft, int n_hop,
                                const Vectorf& window, int original_length = -1) {
    if (X.rows() == 0) {
        return Vectorf();
    }
    
    if (window.size() != n_fft) {
        throw std::invalid_argument("Window size mismatch");
    }
    
    float window_power = window.squaredNorm();
    if (window_power < EPSILON) {
        throw std::runtime_error("Window power is too small");
    }
    
    int64_t n_frames = X.rows();
    int64_t n_f = X.cols();
    int64_t n_fft64 = n_fft;
    int64_t n_hop64 = n_hop;
    
    int64_t expected_length = (n_frames - 1) * n_hop64 + n_fft64;
    
    if (expected_length > 1000000000) {
        throw std::runtime_error("Expected signal length too large");
    }
    
    Matrixcf X_full(n_frames, n_fft);
    X_full.setZero();
    
    X_full.leftCols(n_f) = X;
    
    int nyquist_idx = n_fft / 2;
    bool has_nyquist = (n_fft % 2 == 0);
    
    for (int i = 0; i < n_frames; ++i) {
        if (n_f > 0) {
            X_full(i, 0) = std::real(X_full(i, 0));
        }
        
        if (has_nyquist && n_f > nyquist_idx) {
            X_full(i, nyquist_idx) = std::real(X_full(i, nyquist_idx));
        }
        
        int max_j = has_nyquist ? nyquist_idx : n_f;
        for (int j = 1; j < max_j; ++j) {
            int mirror_idx = n_fft - j;
            if (mirror_idx >= 0 && mirror_idx < n_fft && mirror_idx != j) {
                X_full(i, mirror_idx) = std::conj(X_full(i, j));
            }
        }
    }
    
    Vectorf y = Vectorf::Zero(expected_length);
    Vectorf window_sum = Vectorf::Zero(expected_length);
    
    Eigen::FFT<float> fft;
    
    Vectorcf x_frame_cf(n_fft);
    Vectorf x_frame(n_fft);
    
    for (int i = 0; i < n_frames; ++i) {
        x_frame_cf = X_full.row(i);
        x_frame = fft.inv(x_frame_cf).real();
        
        int64_t start = i * n_hop64;
        
        for (int j = 0; j < n_fft; ++j) {
            if (start + j < expected_length) {
                float windowed_sample = x_frame(j) * window(j);
                y(static_cast<Eigen::Index>(start + j)) += windowed_sample;
                window_sum(static_cast<Eigen::Index>(start + j)) += window(j) * window(j);
            }
        }
    }
    
    for (int i = 0; i < expected_length; ++i) {
        if (window_sum(i) > WINDOW_SUM_EPSILON) {
            y(i) /= window_sum(i);
        }
    }
    
    if (original_length > 0 && original_length < y.size()) {
        y.conservativeResize(original_length);
    }
    
    return y;
}

static Matrixf spectrogram(const Matrixcf& X, float power = 1.0f) {
    if (power <= 0) {
        throw std::invalid_argument("Power must be positive");
    }
    
    Matrixf mag = X.cwiseAbs();
    if (power == 1.0f) {
        return mag;
    } else {
        return mag.array().pow(power);
    }
}

static Matrixf melfilter(int sr, int n_fft, int n_mels, int fmin, int fmax) {
    validate_mel_params(sr, n_fft, n_mels, fmin, fmax);
    
    int n_f = n_fft/2+1;
    Vectorf fft_freqs = (Vectorf::LinSpaced(n_f, 0.f, static_cast<float>(n_f-1))*sr)/n_fft;

    float f_min = 0.f;
    float f_sp = 200.f/3.f;
    float min_log_hz = 1000.f;
    float min_log_mel = (min_log_hz-f_min)/f_sp;
    float logstep = logf(6.4f)/27.f;

    auto hz_to_mel = [=](int hz, bool htk = false) -> float {
        if (htk){
            return 2595.0f*log10f(1.0f+hz/700.0f);
        }
        float mel = (hz-f_min)/f_sp;
        if (hz >= min_log_hz){
            mel = min_log_mel+logf(hz/min_log_hz)/logstep;
        }
        return mel;
    };
    
    auto mel_to_hz = [=](Vectorf &mels, bool htk = false) -> Vectorf {
        if (htk){
            return 700.0f*(Vectorf::Constant(n_mels+2, 10.0f).array().pow(mels.array()/2595.0f)-1.0f);
        }
        return (mels.array()>min_log_mel).select(((mels.array()-min_log_mel)*logstep).exp()*min_log_hz, (mels*f_sp).array()+f_min);
    };

    float min_mel = hz_to_mel(fmin);
    float max_mel = hz_to_mel(fmax);
    Vectorf mels = Vectorf::LinSpaced(n_mels+2, min_mel, max_mel);
    Vectorf mel_f = mel_to_hz(mels);
    Vectorf fdiff = mel_f.segment(1, mel_f.size() - 1) - mel_f.segment(0, mel_f.size() - 1);
    
    for (int i = 0; i < fdiff.size(); ++i) {
        if (std::abs(fdiff(i)) < EPSILON) {
            fdiff(i) = EPSILON;
        }
    }
    
    Matrixf ramps = mel_f.replicate(n_f, 1).transpose().array() - fft_freqs.replicate(n_mels + 2, 1).array();

    Matrixf lower = -ramps.topRows(n_mels).array()/fdiff.segment(0, n_mels).transpose().replicate(1, n_f).array();
    Matrixf upper = ramps.bottomRows(n_mels).array()/fdiff.segment(1, n_mels).transpose().replicate(1, n_f).array();
    Matrixf weights = (lower.array()<upper.array()).select(lower, upper).cwiseMax(0);

    auto enorm = (2.0/(mel_f.segment(2, n_mels)-mel_f.segment(0, n_mels)).array()).transpose().replicate(1, n_f);
    weights = weights.array()*enorm;

    return weights;
}

static Matrixf melspectrogram(const Vectorf& x, int sr, int n_fft, int n_hop,
                             const Vectorf& window, bool center,
                             const std::string& mode, float power,
                             int n_mels, int fmin, int fmax) {
    validate_mel_params(sr, n_fft, n_mels, fmin, fmax);
    
    if (x.empty()) {
        return Matrixf();
    }
    
    Matrixcf X = stft_impl(x, n_fft, n_hop, window, center, mode);
    
    if (X.rows() == 0) {
        throw std::runtime_error("Input signal too short for STFT with given parameters");
    }
    
    Matrixf mel_basis = melfilter(sr, n_fft, n_mels, fmin, fmax);
    
    Matrixf sp = spectrogram(X, power);
    
    Matrixf mel = mel_basis * sp.transpose();
    return mel.transpose();
}

static Matrixf power2db(const Matrixf& x) {
    auto log_sp = 10.0f * x.array().max(EPSILON).log10();
    float max_val = log_sp.maxCoeff();
    return log_sp.cwiseMax(max_val - DB_LIMIT);
}

static Matrixf dct(const Matrixf& x, bool norm, int type) {
    if (x.cols() == 0) {
        return Matrixf();
    }
    
    int N = x.cols();
    Matrixf xi = Matrixf::Zero(N, N);
    xi.rowwise() += Vectorf::LinSpaced(N, 0.f, static_cast<float>(N-1));
    
    Matrixf coeff;
    if (type == 2) {
        coeff = (M_PI * xi.transpose().array() / N * (xi.array() + 0.5)).cos();
        if (!norm) {
            coeff.row(0) *= 0.5f;
        }
    } else if (type == 3) {
        coeff = (M_PI * xi.transpose().array() * (xi.array() + 0.5) / N).cos();
        if (!norm) {
            coeff.row(0) *= 0.5f;
        }
    } else {
        throw std::invalid_argument("Only DCT type 2 and 3 are supported");
    }
    
    Matrixf dct_result = x * coeff.transpose();
    
    if (norm) {
        Vectorf ortho = Vectorf::Constant(N, std::sqrtf(2.0f / N));
        ortho(0) = std::sqrtf(1.0f / N);
        
        for (int i = 0; i < N; ++i) {
            if (std::abs(ortho(i)) < EPSILON) {
                throw std::runtime_error("Division by zero in DCT normalization");
            }
        }
        
        dct_result = dct_result * ortho.asDiagonal();
    }
    
    return dct_result;
}

static Vectorf griffin_lim(const Matrixf& spectrogram, int n_fft, int n_hop,
                           int n_iter = 32, float power = 2.0f,
                           int original_length = -1,
                           bool use_random_init = true,
                           unsigned int seed = 42) {
    
    int n_frames = spectrogram.rows();
    if (n_frames == 0) {
        return Vectorf();
    }
    
    int n_f = spectrogram.cols();
    if (n_f != n_fft / 2 + 1) {
        throw std::invalid_argument("Spectrogram columns must equal n_fft/2 + 1");
    }
    
    Vectorf window = hann_window(n_fft);
    
    auto X = create_matrix<Matrixcf>(n_frames, n_f);
    
    float power_factor = 1.0f / power;
    
    if (use_random_init) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<float> dist(0.0f, 2.0f * M_PI);
        
        for (int i = 0; i < n_frames; ++i) {
            for (int j = 0; j < n_f; ++j) {
                float magnitude = std::pow(std::max(spectrogram(i, j), EPSILON), power_factor);
                float phase = dist(rng);
                (*X)(i, j) = std::polar(magnitude, phase);
            }
        }
    } else {
        for (int i = 0; i < n_frames; ++i) {
            for (int j = 0; j < n_f; ++j) {
                float magnitude = std::pow(std::max(spectrogram(i, j), EPSILON), power_factor);
                (*X)(i, j) = magnitude;
            }
        }
    }
    
    /*
     * Griffin-Lim algorithm implementation:
     * On each iteration:
     * 1. ISTFT to obtain time-domain signal from current spectral estimate X
     * 2. STFT to obtain estimated spectrum X_est
     * 3. Fix magnitudes from target spectrogram while preserving phases from X_est
     * 
     * This is the classic implementation where target magnitudes are enforced
     * on EVERY iteration via std::polar(target_mag, arg(X_est))
     */
    for (int iter = 0; iter < n_iter; ++iter) {
        Vectorf y = istft_with_phases(*X, n_fft, n_hop, window, original_length);
        
        Matrixcf X_est = stft_impl(y, n_fft, n_hop, window, true, "reflect");
        
        int actual_frames = std::min((*X).rows(), X_est.rows());
        int actual_freqs = std::min((*X).cols(), X_est.cols());
        
        for (int i = 0; i < actual_frames; ++i) {
            for (int j = 0; j < actual_freqs; ++j) {
                float target_mag = std::pow(std::max(spectrogram(i, j), EPSILON), power_factor);
                float current_mag = std::abs(X_est(i, j));
                
                if (current_mag > EPSILON) {
                    (*X)(i, j) = std::polar(target_mag, std::arg(X_est(i, j)));
                } else {
                    (*X)(i, j) = target_mag;
                }
            }
        }
    }
    
    return istft_with_phases(*X, n_fft, n_hop, window, original_length);
}

} // namespace internal

class Feature {
public:
    static internal::Matrixcf stft(const std::vector<float>& x, int n_fft, int n_hop,
                                   const std::string& win = "hann",
                                   bool center = true,
                                   const std::string& mode = "reflect") {
        if (x.empty()) {
            return internal::Matrixcf();
        }
        
        internal::Vectorf eigen_x = internal::to_eigen_vector(x);
        internal::Vectorf window = internal::generate_window(win, n_fft);
        
        return internal::stft_impl(eigen_x, n_fft, n_hop, window, center, mode);
    }
    
    static std::vector<float> istft(const internal::Matrixcf& X, int n_fft, int n_hop,
                                    const std::string& win = "hann",
                                    int original_length = -1) {
        if (X.rows() == 0) {
            return std::vector<float>();
        }
        
        internal::Vectorf window = internal::generate_window(win, n_fft);
        
        internal::Vectorf y = internal::istft_with_phases(X, n_fft, n_hop, window, original_length);
        return internal::from_eigen_vector(y);
    }
    
    static internal::Matrixf melspectrogram(const std::vector<float>& x, int sr,
                                           int n_fft, int n_hop,
                                           const std::string& win = "hann",
                                           bool center = true,
                                           const std::string& mode = "reflect",
                                           float power = 2.0f,
                                           int n_mels = 128,
                                           int fmin = 0,
                                           int fmax = -1) {
        if (x.empty()) {
            return internal::Matrixf();
        }
        
        internal::Vectorf eigen_x = internal::to_eigen_vector(x);
        
        if (fmax < 0) {
            fmax = sr / 2;
        }
        
        internal::Vectorf window = internal::generate_window(win, n_fft);
        
        return internal::melspectrogram(eigen_x, sr, n_fft, n_hop, window, center, mode,
                                       power, n_mels, fmin, fmax);
    }
    
    static internal::Matrixf mfcc(const std::vector<float>& x, int sr,
                                 int n_fft, int n_hop,
                                 const std::string& win = "hann",
                                 bool center = true,
                                 const std::string& mode = "reflect",
                                 float power = 2.0f,
                                 int n_mels = 128,
                                 int fmin = 0,
                                 int fmax = -1,
                                 int n_mfcc = 20,
                                 bool norm = true,
                                 int type = 2) {
        if (x.empty()) {
            return internal::Matrixf();
        }
        
        if (n_mfcc > n_mels) {
            throw std::invalid_argument("n_mfcc cannot exceed n_mels");
        }
        
        if (type != 2 && type != 3) {
            throw std::invalid_argument("DCT type must be 2 or 3");
        }
        
        internal::Vectorf eigen_x = internal::to_eigen_vector(x);
        
        if (fmax < 0) {
            fmax = sr / 2;
        }
        
        internal::Vectorf window = internal::generate_window(win, n_fft);
        
        internal::Matrixf mel = internal::melspectrogram(eigen_x, sr, n_fft, n_hop,
                                                        window, center, mode,
                                                        power, n_mels, fmin, fmax);
        
        if (mel.rows() == 0) {
            return internal::Matrixf(0, n_mfcc);
        }
        
        internal::Matrixf mel_db = internal::power2db(mel);
        
        internal::Matrixf dct_result = internal::dct(mel_db, norm, type);
        
        return dct_result.leftCols(n_mfcc);
    }
};

class GriffinLim {
public:
    static std::vector<float> reconstruct(
        const internal::Matrixf& spectrogram,
        int n_fft,
        int n_hop,
        int n_iter = 32,
        float power = 2.0f,
        int original_length = -1,
        bool use_random_init = true,
        unsigned int seed = 42) {
        
        if (spectrogram.rows() == 0) {
            return std::vector<float>();
        }
        
        internal::Vectorf y = internal::griffin_lim(spectrogram, n_fft, n_hop,
                                                   n_iter, power,
                                                   original_length,
                                                   use_random_init, seed);
        
        return internal::from_eigen_vector(y);
    }
    
    static std::vector<float> reconstruct(
        const SpectrogramData& spectrogram,
        int n_fft,
        int n_hop,
        int n_iter = 32,
        float power = 2.0f,
        int original_length = -1,
        bool use_random_init = true,
        unsigned int seed = 42) {
        
        if (spectrogram.n_frames == 0) {
            return std::vector<float>();
        }
        
        internal::Matrixf S = spectrogram.to_eigen_mag();
        internal::Vectorf y = internal::griffin_lim(S, n_fft, n_hop,
                                                   n_iter, power,
                                                   original_length,
                                                   use_random_init, seed);
        
        return internal::from_eigen_vector(y);
    }
};

} // namespace librosa

#endif