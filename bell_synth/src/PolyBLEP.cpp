// =============================================================================
// PolyBLEP.cpp - РЕАЛИЗАЦИЯ ОСЦИЛЛЯТОРА С АНТИАЛИАСИНГОМ
// =============================================================================
#include "PolyBLEP.h"
#include <cmath>
#include <cstdint>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static constexpr double TWO_PI = 2.0 * M_PI;

template<typename T>
inline T square_number(const T &x) {
    return x * x;
}

inline double blep(double t, double dt) {
    if (t < dt) {
        return -square_number(t / dt - 1);
    } else if (t > 1 - dt) {
        return square_number((t - 1) / dt + 1);
    } else {
        return 0;
    }
}

inline double blamp(double t, double dt) {
    if (t < dt) {
        t = t / dt - 1;
        return -1 / 3.0 * square_number(t) * t;
    } else if (t > 1 - dt) {
        t = (t - 1) / dt + 1;
        return 1 / 3.0 * square_number(t) * t;
    } else {
        return 0;
    }
}

template<typename T>
inline int64_t bitwiseOrZero(const T &t) {
    return static_cast<int64_t>(t) | 0;
}

PolyBLEP::PolyBLEP(double sampleRate)
        : sampleRate_(sampleRate)
        , amplitude_(1.0)
        , t_(0.0)
        , waveform_(SINE)
        , pulseWidth_(0.5) {
    setSampleRate(sampleRate);
    setFrequency(440.0);
    setWaveform(SINE);
    setPulseWidth(0.5);
}

void PolyBLEP::setdt(double time) {
    freqInSecondsPerSample_ = time;
}

void PolyBLEP::setFrequency(double freqInHz) {
    setdt(freqInHz / sampleRate_);
}

void PolyBLEP::setSampleRate(double sampleRate) {
    const double freqInHz = getFreqInHz();
    this->sampleRate_ = sampleRate;
    setFrequency(freqInHz);
}

double PolyBLEP::getFreqInHz() const {
    return freqInSecondsPerSample_ * sampleRate_;
}

void PolyBLEP::setPulseWidth(double pulseWidth) {
    pulseWidth_ = std::clamp(pulseWidth, 0.0001, 0.9999);
}

void PolyBLEP::sync(double phase) {
    t_ = phase;
    if (t_ >= 0) {
        t_ -= bitwiseOrZero(t_);
    } else {
        t_ += 1 - bitwiseOrZero(t_);
    }
}

void PolyBLEP::setWaveform(Waveform waveform) {
    waveform_ = waveform;
}

double PolyBLEP::get() const {
    if (getFreqInHz() >= sampleRate_ / 4) {
        return sin();
    }
    
    switch (waveform_) {
        case SINE: return sin();
        case COSINE: return cos();
        case TRIANGLE: return tri();
        case SQUARE: return sqr();
        case RECTANGLE: return rect();
        case SAWTOOTH: return saw();
        case RAMP: return ramp();
        case MODIFIED_TRIANGLE: return tri2();
        case MODIFIED_SQUARE: return sqr2();
        case HALF_WAVE_RECTIFIED_SINE: return half();
        case FULL_WAVE_RECTIFIED_SINE: return full();
        case TRIANGULAR_PULSE: return trip();
        case TRAPEZOID_FIXED: return trap();
        case TRAPEZOID_VARIABLE: return trap2();
        default: return 0.0;
    }
}

void PolyBLEP::inc() {
    t_ += freqInSecondsPerSample_;
    t_ -= bitwiseOrZero(t_);
}

double PolyBLEP::getAndInc() {
    const double sample = get();
    inc();
    return sample;
}

double PolyBLEP::sin() const {
    return amplitude_ * std::sin(TWO_PI * t_);
}

double PolyBLEP::cos() const {
    return amplitude_ * std::cos(TWO_PI * t_);
}

double PolyBLEP::half() const {
    double t2 = t_ + 0.5;
    t2 -= bitwiseOrZero(t2);

    double y = (t_ < 0.5 ? 2 * std::sin(TWO_PI * t_) - 2 / M_PI : -2 / M_PI);
    y += TWO_PI * freqInSecondsPerSample_ * (blamp(t_, freqInSecondsPerSample_) + blamp(t2, freqInSecondsPerSample_));

    return amplitude_ * y;
}

double PolyBLEP::full() const {
    double _t = t_ + 0.25;
    _t -= bitwiseOrZero(_t);

    double y = 2 * std::sin(M_PI * _t) - 4 / M_PI;
    y += TWO_PI * freqInSecondsPerSample_ * blamp(_t, freqInSecondsPerSample_);

    return amplitude_ * y;
}

double PolyBLEP::tri() const {
    double t1 = t_ + 0.25;
    t1 -= bitwiseOrZero(t1);

    double t2 = t_ + 0.75;
    t2 -= bitwiseOrZero(t2);

    double y = t_ * 4;

    if (y >= 3) {
        y -= 4;
    } else if (y > 1) {
        y = 2 - y;
    }

    y += 4 * freqInSecondsPerSample_ * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    return amplitude_ * y;
}

double PolyBLEP::tri2() const {
    double pulseWidth = std::clamp(pulseWidth_, 0.0001, 0.9999);

    double t1 = t_ + 0.5 * pulseWidth;
    t1 -= bitwiseOrZero(t1);

    double t2 = t_ + 1 - 0.5 * pulseWidth;
    t2 -= bitwiseOrZero(t2);

    double y = t_ * 2;

    if (y >= 2 - pulseWidth) {
        y = (y - 2) / pulseWidth;
    } else if (y >= pulseWidth) {
        y = 1 - (y - pulseWidth) / (1 - pulseWidth);
    } else {
        y /= pulseWidth;
    }

    y += freqInSecondsPerSample_ / (pulseWidth - pulseWidth * pulseWidth) * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    return amplitude_ * y;
}

double PolyBLEP::trip() const {
    double t1 = t_ + 0.75 + 0.5 * pulseWidth_;
    t1 -= bitwiseOrZero(t1);

    double y;
    if (t1 >= pulseWidth_) {
        y = -pulseWidth_;
    } else {
        y = 4 * t1;
        y = (y >= 2 * pulseWidth_ ? 4 - y / pulseWidth_ - pulseWidth_ : y / pulseWidth_ - pulseWidth_);
    }

    if (pulseWidth_ > 0) {
        double t2 = t1 + 1 - 0.5 * pulseWidth_;
        t2 -= bitwiseOrZero(t2);

        double t3 = t1 + 1 - pulseWidth_;
        t3 -= bitwiseOrZero(t3);
        y += 2 * freqInSecondsPerSample_ / pulseWidth_ * (blamp(t1, freqInSecondsPerSample_) - 2 * blamp(t2, freqInSecondsPerSample_) + blamp(t3, freqInSecondsPerSample_));
    }
    return amplitude_ * y;
}

double PolyBLEP::trap() const {
    double y = 4 * t_;
    if (y >= 3) {
        y -= 4;
    } else if (y > 1) {
        y = 2 - y;
    }
    y = std::clamp(2 * y, -1.0, 1.0);

    double t1 = t_ + 0.125;
    t1 -= bitwiseOrZero(t1);

    double t2 = t1 + 0.5;
    t2 -= bitwiseOrZero(t2);

    y += 4 * freqInSecondsPerSample_ * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    t1 = t_ + 0.375;
    t1 -= bitwiseOrZero(t1);

    t2 = t1 + 0.5;
    t2 -= bitwiseOrZero(t2);

    y += 4 * freqInSecondsPerSample_ * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    return amplitude_ * y;
}

double PolyBLEP::trap2() const {
    double pulseWidth = std::clamp(pulseWidth_, 0.0001, 0.9999);
    double scale = 1 / (1 - pulseWidth);

    double y = 4 * t_;
    if (y >= 3) {
        y -= 4;
    } else if (y > 1) {
        y = 2 - y;
    }
    y = std::clamp(scale * y, -1.0, 1.0);

    double t1 = t_ + 0.25 - 0.25 * pulseWidth;
    t1 -= bitwiseOrZero(t1);

    double t2 = t1 + 0.5;
    t2 -= bitwiseOrZero(t2);

    y += scale * 2 * freqInSecondsPerSample_ * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    t1 = t_ + 0.25 + 0.25 * pulseWidth;
    t1 -= bitwiseOrZero(t1);

    t2 = t1 + 0.5;
    t2 -= bitwiseOrZero(t2);

    y += scale * 2 * freqInSecondsPerSample_ * (blamp(t1, freqInSecondsPerSample_) - blamp(t2, freqInSecondsPerSample_));

    return amplitude_ * y;
}

double PolyBLEP::sqr() const {
    double t2 = t_ + 0.5;
    t2 -= bitwiseOrZero(t2);

    double y = t_ < 0.5 ? 1 : -1;
    y += blep(t_, freqInSecondsPerSample_) - blep(t2, freqInSecondsPerSample_);

    return amplitude_ * y;
}

double PolyBLEP::sqr2() const {
    double t1 = t_ + 0.875 + 0.25 * (pulseWidth_ - 0.5);
    t1 -= bitwiseOrZero(t1);

    double t2 = t_ + 0.375 + 0.25 * (pulseWidth_ - 0.5);
    t2 -= bitwiseOrZero(t2);

    double y = t1 < 0.5 ? 1 : -1;
    y += blep(t1, freqInSecondsPerSample_) - blep(t2, freqInSecondsPerSample_);

    t1 += 0.5 * (1 - pulseWidth_);
    t1 -= bitwiseOrZero(t1);

    t2 += 0.5 * (1 - pulseWidth_);
    t2 -= bitwiseOrZero(t2);

    y += t1 < 0.5 ? 1 : -1;
    y += blep(t1, freqInSecondsPerSample_) - blep(t2, freqInSecondsPerSample_);

    return amplitude_ * 0.5 * y;
}

double PolyBLEP::rect() const {
    double t2 = t_ + 1 - pulseWidth_;
    t2 -= bitwiseOrZero(t2);

    double y = -2 * pulseWidth_;
    if (t_ < pulseWidth_) {
        y += 2;
    }

    y += blep(t_, freqInSecondsPerSample_) - blep(t2, freqInSecondsPerSample_);

    return amplitude_ * y;
}

double PolyBLEP::saw() const {
    double _t = t_ + 0.5;
    _t -= bitwiseOrZero(_t);

    double y = 2 * _t - 1;
    y -= blep(_t, freqInSecondsPerSample_);

    return amplitude_ * y;
}

double PolyBLEP::ramp() const {
    double _t = t_;
    _t -= bitwiseOrZero(_t);

    double y = 1 - 2 * _t;
    y += blep(_t, freqInSecondsPerSample_);

    return amplitude_ * y;
}