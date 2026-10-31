#include "cylindrical_strain.hpp"
#include <cmath>
#include <stdexcept>

namespace BellOptimizer::Math {
CylindricalStrain cylindrical_small_strain(const CylindricalDisplacementJet& u) {
    if (!(u.r > 0.0) || !std::isfinite(u.r))
        throw std::invalid_argument("cylindrical strain: radius must be finite and positive");
    const double invr = 1.0 / u.r;
    return {{
        u.dur_dr,
        invr * u.duphi_dphi + invr * u.ur,
        u.duz_dz,
        invr * u.dur_dphi + u.duphi_dr - invr * u.uphi,
        u.dur_dz + u.duz_dr,
        u.duphi_dz + invr * u.duz_dphi
    }};
}

ElasticityMatrix6 isotropic_elasticity_3d(double E, double nu) {
    if (!(E > 0.0) || !std::isfinite(E) || !(nu > -1.0 && nu < 0.5) || !std::isfinite(nu))
        throw std::invalid_argument("isotropic elasticity: require E>0 and -1<nu<0.5");
    ElasticityMatrix6 c{};
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = lambda;
        c[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] += 2.0 * mu;
    }
    c[3][3] = c[4][4] = c[5][5] = mu;
    return c;
}

double elastic_strain_energy_density(const CylindricalStrain& e, const ElasticityMatrix6& c) {
    double energy2 = 0.0;
    for (std::size_t i = 0; i < e.size(); ++i)
        for (std::size_t j = 0; j < e.size(); ++j)
            energy2 += e[i] * c[i][j] * e[j];
    if (!std::isfinite(energy2) || energy2 < -1e-10)
        throw std::runtime_error("elastic energy density: invalid constitutive result");
    return 0.5 * (energy2 < 0.0 ? 0.0 : energy2);
}
} // namespace BellOptimizer::Math
