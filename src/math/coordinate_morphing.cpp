#include "coordinate_morphing.hpp"

#include <limits>

namespace BellOptimizer::Math {

AcousticWaveVector3D transform_wave_vector(const AcousticWaveVector3D& k,
                                           const CoordinateScale3D& scale) {
    scale.validate();
    return {k.kx / scale.x, k.ky / scale.y, k.kz / scale.z};
}

double wave_number(const AcousticWaveVector3D& k) {
    return std::sqrt(k.kx * k.kx + k.ky * k.ky + k.kz * k.kz);
}

double transformed_wave_number(const AcousticWaveVector3D& k,
                               const CoordinateScale3D& scale) {
    return wave_number(transform_wave_vector(k, scale));
}

double transformed_acoustic_frequency_hz(const AcousticWaveVector3D& k,
                                         const CoordinateScale3D& scale,
                                         double sound_speed_m_s) {
    if (!(sound_speed_m_s > 0.0) || !std::isfinite(sound_speed_m_s))
        throw std::invalid_argument("sound speed must be finite and positive");
    const double k_prime = transformed_wave_number(k, scale);
    return sound_speed_m_s * k_prime / (2.0 * 3.14159265358979323846);
}

} // namespace BellOptimizer::Math
