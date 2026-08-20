// =============================================================================
// PolyBLEP.h - ОСЦИЛЛЯТОР С АНТИАЛИАСИНГОМ
// =============================================================================
#ifndef POLY_BLEP_H_INCLUDED
#define POLY_BLEP_H_INCLUDED

#include <cmath>
#include <cstdint>
#include <algorithm>

class PolyBLEP {
public:
    enum Waveform {
        SINE = 0,
        COSINE,
        TRIANGLE,
        SQUARE,
        RECTANGLE,
        SAWTOOTH,
        RAMP,
        MODIFIED_TRIANGLE,
        MODIFIED_SQUARE,
        HALF_WAVE_RECTIFIED_SINE,
        FULL_WAVE_RECTIFIED_SINE,
        TRIANGULAR_PULSE,
        TRAPEZOID_FIXED,
        TRAPEZOID_VARIABLE
    };

    PolyBLEP(double sampleRate = 44100.0);
    ~PolyBLEP() = default;

    void setFrequency(double freqInHz);
    void setSampleRate(double sampleRate);
    void setWaveform(Waveform waveform);
    void setPulseWidth(double pw);
    void setAmplitude(double amp) { amplitude_ = amp; }
    double getAmplitude() const { return amplitude_; }

    double get() const;
    void inc();
    double getAndInc();
    double getFreqInHz() const;
    void sync(double phase);

protected:
    Waveform waveform_;
    double sampleRate_;
    double freqInSecondsPerSample_;
    double amplitude_;
    double pulseWidth_;
    double t_;

    void setdt(double time);
    double sin() const;
    double cos() const;
    double half() const;
    double full() const;
    double tri() const;
    double tri2() const;
    double trip() const;
    double trap() const;
    double trap2() const;
    double sqr() const;
    double sqr2() const;
    double rect() const;
    double saw() const;
    double ramp() const;
};

#endif // POLY_BLEP_H_INCLUDED