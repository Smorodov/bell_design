#include "directional_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Acoustics {
namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kSoundSpeed = 343.0;

bool validSurface(const ModalMode& mode, const std::vector<Point2D>& profile) {
    const std::size_t expected = profile.size() * mode.outer_surface_harmonics.size();
    return profile.size() >= 2 && !mode.outer_surface_harmonics.empty() &&
        mode.outer_surface_fourier_radial.size() == expected &&
        mode.outer_surface_fourier_axial.size() == expected &&
        std::isfinite(mode.frequency_hz) && mode.frequency_hz > 0.0;
}

double basisValue(const std::vector<double>& values, std::size_t point,
                  std::size_t basisCount, std::size_t basis) {
    return values[point * basisCount + basis];
}

double fourierValue(const std::vector<double>& values, std::size_t point,
                    const std::vector<int>& labels, double angle) {
    double result = 0.0;
    for (std::size_t b = 0; b < labels.size(); ++b) {
        const int label = labels[b];
        const double coefficient = values[point * labels.size() + b];
        if (label == 0) result += coefficient;
        else if (label > 0) result += coefficient * std::cos(label * angle);
        else result += coefficient * std::sin(-label * angle);
    }
    return result;
}

} // namespace

std::complex<double> evaluateDirectionalRadiationDirect(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double listenerPolarRad, double listenerAzimuthRad, unsigned azimuthSamples) {
    if (!validSurface(mode, profile) || !std::isfinite(listenerPolarRad) ||
        !std::isfinite(listenerAzimuthRad) || azimuthSamples < 4)
        throw std::invalid_argument("invalid input for direct directional radiation");

    const double waveNumber = 2.0 * kPi * mode.frequency_hz / kSoundSpeed;
    const double sinEl = std::sin(listenerPolarRad);
    const double cosEl = std::cos(listenerPolarRad);
    std::complex<double> integral{};
    for (std::size_t j = 0; j + 1 < profile.size(); ++j) {
        const double dr = profile[j + 1].r - profile[j].r;
        const double dz = profile[j + 1].z - profile[j].z;
        const double ds = std::hypot(dr, dz);
        if (!(ds > 0.0)) continue;
        const double nr = dz / ds, nz = -dr / ds;
        const double r = 0.5 * (profile[j].r + profile[j + 1].r) * 1.0e-3;
        const double z = 0.5 * (profile[j].z + profile[j + 1].z) * 1.0e-3;
        const double dsMeters = ds * 1.0e-3;
        for (unsigned a = 0; a < azimuthSamples; ++a) {
            const double phi = 2.0 * kPi * (static_cast<double>(a) + 0.5) /
                               static_cast<double>(azimuthSamples);
            const double ur = 0.5 * (
                fourierValue(mode.outer_surface_fourier_radial, j,
                             mode.outer_surface_harmonics, phi) +
                fourierValue(mode.outer_surface_fourier_radial, j + 1,
                             mode.outer_surface_harmonics, phi));
            const double uz = 0.5 * (
                fourierValue(mode.outer_surface_fourier_axial, j,
                             mode.outer_surface_harmonics, phi) +
                fourierValue(mode.outer_surface_fourier_axial, j + 1,
                             mode.outer_surface_harmonics, phi));
            const double normal = ur * nr + uz * nz;
            const double phase = -waveNumber *
                (r * sinEl * std::cos(phi - listenerAzimuthRad) + z * cosEl);
            integral += normal * r * dsMeters * (2.0 * kPi / azimuthSamples) *
                        std::complex<double>(std::cos(phase), std::sin(phase));
        }
    }
    return integral;
}

AcousticDirectionalField buildDirectionalRadiationField(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double listenerPolarRad, double observationRadiusM) {
    if (!validSurface(mode, profile) || !std::isfinite(listenerPolarRad) ||
        !(std::isfinite(observationRadiusM) && observationRadiusM > 0.0))
        throw std::invalid_argument("invalid modal surface field for directional radiation");

    const double waveNumber = 2.0 * kPi * mode.frequency_hz / kSoundSpeed;
    unsigned maxOrder = 0;
    for (int label : mode.outer_surface_harmonics)
        maxOrder = std::max(maxOrder, static_cast<unsigned>(std::abs(label)));
    std::vector<std::complex<double>> cosine(maxOrder + 1), sine(maxOrder + 1);
    const double sinEl = std::sin(listenerPolarRad), cosEl = std::cos(listenerPolarRad);
    const std::size_t basisCount = mode.outer_surface_harmonics.size();

    for (std::size_t j = 0; j + 1 < profile.size(); ++j) {
        const double dr = profile[j + 1].r - profile[j].r;
        const double dz = profile[j + 1].z - profile[j].z;
        const double ds = std::hypot(dr, dz);
        if (!(ds > 0.0)) continue;
        const double nr = dz / ds, nz = -dr / ds;
        const double r = 0.5 * (profile[j].r + profile[j + 1].r) * 1.0e-3;
        const double z = 0.5 * (profile[j].z + profile[j + 1].z) * 1.0e-3;
        const double dsMeters = ds * 1.0e-3;
        const double x = waveNumber * r * sinEl;
        const std::complex<double> axialPhase = std::polar(1.0, -waveNumber * z * cosEl);
        const double surfaceWeight = r * dsMeters;
        for (std::size_t b = 0; b < basisCount; ++b) {
            const int label = mode.outer_surface_harmonics[b];
            const unsigned order = static_cast<unsigned>(std::abs(label));
            const double ur = 0.5 * (basisValue(mode.outer_surface_fourier_radial, j, basisCount, b) +
                                      basisValue(mode.outer_surface_fourier_radial, j + 1, basisCount, b));
            const double uz = 0.5 * (basisValue(mode.outer_surface_fourier_axial, j, basisCount, b) +
                                      basisValue(mode.outer_surface_fourier_axial, j + 1, basisCount, b));
            const double normal = ur * nr + uz * nz;
            const double bessel = std::cyl_bessel_j(static_cast<double>(order), x);
            const auto angularFactor = 2.0 * kPi *
                std::pow(std::complex<double>(0.0, -1.0), static_cast<int>(order)) * bessel;
            const auto contribution = normal * surfaceWeight * axialPhase * angularFactor;
            if (label == 0) cosine[0] += contribution;
            else if (label > 0) cosine[order] += contribution;
            else sine[order] += contribution;
        }
    }

    AcousticDirectionalField field;
    field.observation_radius_m = observationRadiusM;
    field.observation_polar_rad = listenerPolarRad;
    field.orders.reserve(2 * maxOrder + 1);
    field.coefficients.reserve(2 * maxOrder + 1);
    for (int order = -static_cast<int>(maxOrder); order <= static_cast<int>(maxOrder); ++order) {
        field.orders.push_back(order);
        if (order == 0) field.coefficients.push_back(cosine[0]);
        else {
            const auto m = static_cast<unsigned>(std::abs(order));
            field.coefficients.push_back(order > 0
                ? 0.5 * (cosine[m] - std::complex<double>(0.0, 1.0) * sine[m])
                : 0.5 * (cosine[m] + std::complex<double>(0.0, 1.0) * sine[m]));
        }
    }
    if (!field.valid()) throw std::runtime_error("constructed invalid directional radiation field");
    return field;
}


bool AcousticPolarDirectionalField::valid() const {
    if (polar_angles_rad.empty() || polar_angles_rad.size() != fields.size()) return false;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (!fields[i].valid() || !std::isfinite(polar_angles_rad[i]) ||
            polar_angles_rad[i] < 0.0 || polar_angles_rad[i] > kPi) return false;
        if (i > 0 && !(polar_angles_rad[i] > polar_angles_rad[i - 1])) return false;
        if (std::abs(fields[i].observation_polar_rad - polar_angles_rad[i]) > 1.0e-12) return false;
        if (i > 0 && (fields[i].orders != fields[0].orders ||
                      fields[i].observation_radius_m != fields[0].observation_radius_m)) return false;
    }
    return true;
}

namespace {
std::complex<double> evaluatePolarField(const AcousticPolarDirectionalField& polarField,
                                        double polarRad, double azimuthRad,
                                        double frameRotationRad) {
    if (!polarField.valid() || !std::isfinite(polarRad) || polarRad < 0.0 ||
        polarRad > kPi || !std::isfinite(azimuthRad) || !std::isfinite(frameRotationRad))
        throw std::invalid_argument("invalid polar directional field query");
    if (polarField.fields.size() == 1) {
        if (std::abs(polarRad - polarField.polar_angles_rad.front()) > 1.0e-10)
            throw std::out_of_range("polar query lies outside single-node directional field");
        return polarField.fields.front().evaluateRotated(azimuthRad, frameRotationRad);
    }
    if (polarRad < polarField.polar_angles_rad.front() - 1.0e-12 ||
        polarRad > polarField.polar_angles_rad.back() + 1.0e-12)
        throw std::out_of_range("polar query lies outside precomputed directional field range");
    const auto upper = std::lower_bound(polarField.polar_angles_rad.begin(),
                                        polarField.polar_angles_rad.end(), polarRad);
    if (upper == polarField.polar_angles_rad.begin())
        return polarField.fields.front().evaluateRotated(azimuthRad, frameRotationRad);
    if (upper == polarField.polar_angles_rad.end())
        return polarField.fields.back().evaluateRotated(azimuthRad, frameRotationRad);
    const std::size_t hi = static_cast<std::size_t>(upper - polarField.polar_angles_rad.begin());
    const std::size_t lo = hi - 1;
    const double t = (polarRad - polarField.polar_angles_rad[lo]) /
                     (polarField.polar_angles_rad[hi] - polarField.polar_angles_rad[lo]);
    const auto a = polarField.fields[lo].evaluateRotated(azimuthRad, frameRotationRad);
    const auto b = polarField.fields[hi].evaluateRotated(azimuthRad, frameRotationRad);
    return (1.0 - t) * a + t * b;
}
} // namespace

std::complex<double> AcousticPolarDirectionalField::evaluate(double polarRad,
                                                              double azimuthRad) const {
    return evaluatePolarField(*this, polarRad, azimuthRad, 0.0);
}

std::complex<double> AcousticPolarDirectionalField::evaluateRotated(
    double polarRad, double azimuthRad, double frameRotationRad) const {
    return evaluatePolarField(*this, polarRad, azimuthRad, frameRotationRad);
}

AcousticPolarDirectionalField buildDirectionalRadiationPolarField(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double observationRadiusM, unsigned polarSamples) {
    if (polarSamples < 3 || polarSamples > 721)
        throw std::invalid_argument("polar radiation grid must contain 3..721 samples");
    AcousticPolarDirectionalField result;
    result.polar_angles_rad.reserve(polarSamples);
    result.fields.reserve(polarSamples);
    for (unsigned i = 0; i < polarSamples; ++i) {
        const double polar = kPi * static_cast<double>(i) /
                            static_cast<double>(polarSamples - 1);
        result.polar_angles_rad.push_back(polar);
        result.fields.push_back(buildDirectionalRadiationField(mode, profile, polar,
                                                                observationRadiusM));
    }
    if (!result.valid()) throw std::runtime_error("constructed invalid polar directional radiation field");
    return result;
}


AcousticPolarDirectionalField buildDirectionalRadiationPolarFieldRange(
    const ModalMode& mode, const std::vector<Point2D>& profile,
    double observationRadiusM, double minimumPolarRad,
    double maximumPolarRad, unsigned polarSamples) {
    if (!std::isfinite(minimumPolarRad) || !std::isfinite(maximumPolarRad) ||
        minimumPolarRad < 0.0 || maximumPolarRad > kPi ||
        maximumPolarRad < minimumPolarRad || polarSamples < 1 || polarSamples > 721 ||
        (maximumPolarRad > minimumPolarRad && polarSamples < 2) ||
        (maximumPolarRad == minimumPolarRad && polarSamples != 1))
        throw std::invalid_argument("invalid polar interval or sample count");
    AcousticPolarDirectionalField result;
    result.polar_angles_rad.reserve(polarSamples);
    result.fields.reserve(polarSamples);
    for (unsigned i = 0; i < polarSamples; ++i) {
        const double fraction = polarSamples == 1 ? 0.0 :
            static_cast<double>(i) / static_cast<double>(polarSamples - 1);
        const double polar = minimumPolarRad + fraction * (maximumPolarRad - minimumPolarRad);
        result.polar_angles_rad.push_back(polar);
        result.fields.push_back(buildDirectionalRadiationField(mode, profile, polar,
                                                                observationRadiusM));
    }
    if (!result.valid()) throw std::runtime_error("constructed invalid local polar directional field");
    return result;
}



std::vector<StereoRadiationTransfer> PreparedStereoRadiation::query(
    double headYawRad, double centreAzimuthOffsetRad) const {
    if (fields_.empty() || !std::isfinite(headYawRad) || !std::isfinite(centreAzimuthOffsetRad) ||
        std::abs(headYawRad) > maximumAbsYawRad_ + 1.0e-12)
        throw std::out_of_range("head yaw is outside the prepared stereo radiation envelope");
    const auto listener = makeStereoListenerGeometry(
        centreDistanceM_, centrePolarRad_, centreAzimuthRad_ + centreAzimuthOffsetRad,
        earSpacingM_, headYawRad);
    std::vector<StereoRadiationTransfer> result;
    result.reserve(fields_.size());
    for (const auto& field : fields_) {
        StereoRadiationTransfer transfer;
        transfer.left = field.evaluate(listener.left.polar_rad, listener.left.azimuth_rad);
        transfer.right = field.evaluate(listener.right.polar_rad, listener.right.azimuth_rad);
        transfer.left_range_scale = centreDistanceM_ / listener.left.distance_m;
        transfer.right_range_scale = centreDistanceM_ / listener.right.distance_m;
        transfer.left_delay_s = listener.left.propagation_delay_s;
        transfer.right_delay_s = listener.right.propagation_delay_s;
        result.push_back(transfer);
    }
    return result;
}

PreparedStereoRadiation prepareStereoRadiation(
    const std::vector<const ModalMode*>& modes,
    const std::vector<Point2D>& profile,
    double centreDistanceM, double centrePolarRad,
    double centreAzimuthRad, double earSpacingM,
    double maximumAbsYawRad, double maximumPolarStepRad) {
    if (modes.empty() || profile.size() < 2 ||
        !(std::isfinite(centreDistanceM) && centreDistanceM > 0.0) ||
        !std::isfinite(centrePolarRad) || centrePolarRad < 0.0 || centrePolarRad > kPi ||
        !std::isfinite(centreAzimuthRad) ||
        !(std::isfinite(earSpacingM) && earSpacingM >= 0.0) ||
        !(std::isfinite(maximumAbsYawRad) && maximumAbsYawRad >= 0.0 && maximumAbsYawRad <= 0.5*kPi) ||
        !(std::isfinite(maximumPolarStepRad) && maximumPolarStepRad > 0.0 && maximumPolarStepRad <= kPi))
        throw std::invalid_argument("invalid prepared stereo radiation configuration");

    const auto yawMinus = makeStereoListenerGeometry(centreDistanceM, centrePolarRad,
        centreAzimuthRad, earSpacingM, -maximumAbsYawRad);
    const auto yawPlus = makeStereoListenerGeometry(centreDistanceM, centrePolarRad,
        centreAzimuthRad, earSpacingM, maximumAbsYawRad);
    const double polarMin = std::min({yawMinus.left.polar_rad, yawMinus.right.polar_rad,
                                      yawPlus.left.polar_rad, yawPlus.right.polar_rad});
    const double polarMax = std::max({yawMinus.left.polar_rad, yawMinus.right.polar_rad,
                                      yawPlus.left.polar_rad, yawPlus.right.polar_rad});
    const double span = polarMax - polarMin;
    const unsigned nodes = span < 1.0e-12 ? 1U :
        std::max(2U, static_cast<unsigned>(std::ceil(span / maximumPolarStepRad)) + 1U);

    PreparedStereoRadiation prepared;
    prepared.centreDistanceM_ = centreDistanceM;
    prepared.centrePolarRad_ = centrePolarRad;
    prepared.centreAzimuthRad_ = centreAzimuthRad;
    prepared.earSpacingM_ = earSpacingM;
    prepared.maximumAbsYawRad_ = maximumAbsYawRad;
    prepared.fields_.reserve(modes.size());
    prepared.mode_frequencies_hz_.reserve(modes.size());
    for (const auto* mode : modes) {
        if (!mode) throw std::invalid_argument("null mode in prepared stereo radiation input");
        prepared.fields_.push_back(buildDirectionalRadiationPolarFieldRange(
            *mode, profile, centreDistanceM, polarMin, polarMax, nodes));
        prepared.mode_frequencies_hz_.push_back(mode->frequency_hz);
        prepared.surface_integral_build_count_ += nodes;
    }
    return prepared;
}

} // namespace BellOptimizer::Acoustics
