#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Math {

// Positive diagonal spatial scale. A value below one compresses that axis;
// a value above one stretches it. This layer describes a coordinate morph,
// not a new FEM solve.
struct CoordinateScale3D {
    double x = 1.0;
    double y = 1.0;
    double z = 1.0;

    void validate() const {
        if (!(x > 0.0) || !(y > 0.0) || !(z > 0.0) ||
            !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            throw std::invalid_argument("CoordinateScale3D requires finite positive scales");
    }

    double volume_scale() const { validate(); return x * y * z; }

    std::array<double, 3> point(const std::array<double, 3>& p) const {
        validate();
        return {x * p[0], y * p[1], z * p[2]};
    }
};

struct AcousticWaveVector3D {
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
};

// Pulling a scalar acoustic field p'(x') = p(S^{-1}x') transforms the
// wave-vector covariantly: k' = S^{-T} k. For diagonal scaling this is simply
// k_i' = k_i / s_i.
AcousticWaveVector3D transform_wave_vector(const AcousticWaveVector3D& k,
                                           const CoordinateScale3D& scale);

double wave_number(const AcousticWaveVector3D& k);

double transformed_wave_number(const AcousticWaveVector3D& k,
                               const CoordinateScale3D& scale);

// Exact for a plane-wave component in an isotropic, homogeneous acoustic
// medium. It is the first reduced-order frequency correction used by the
// coordinate-morphing layer; it is not a replacement for a full 3-D solve.
double transformed_acoustic_frequency_hz(const AcousticWaveVector3D& k,
                                         const CoordinateScale3D& scale,
                                         double sound_speed_m_s);

} // namespace BellOptimizer::Math
