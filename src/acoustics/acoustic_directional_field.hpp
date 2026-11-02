#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace BellOptimizer::Acoustics {

// Compact immutable-base representation of a directional far-field transfer
// function at one declared observation condition. Coefficients are indexed by
// azimuthal order m and are complex so that propagation/radiation phase is not
// discarded.
//
// This type intentionally does not compute radiation coefficients. It only
// defines the representation and the mathematically exact azimuthal rotation
// operation. Population from FEM/acoustic radiation remains a separate gate.
struct AcousticDirectionalField {
    double observation_radius_m = 1.0;
    double observation_elevation_rad = 0.0;
    std::vector<int> orders;
    std::vector<std::complex<double>> coefficients;

    bool valid() const;
    std::size_t size() const { return orders.size(); }

    // Return a rotated observation of the same physical base field. Positive
    // theta follows the declared positive azimuth convention of the owning
    // acoustic model; the representation itself applies exp(-i*m*theta).
    AcousticDirectionalField rotated(double theta_rad) const;

    // Evaluate the retained harmonic field at azimuth theta.
    std::complex<double> evaluate(double theta_rad) const;
};

} // namespace BellOptimizer::Acoustics
