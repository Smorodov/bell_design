#pragma once

#include <array>

namespace BellOptimizer::Math {

// Cylindrical displacement field and its partial derivatives at (r, phi, z).
// Derivatives are physical-coordinate derivatives, not derivatives with
// respect to a normalized/ computational coordinate.
struct CylindricalDisplacementJet {
    double r = 0.0;
    double ur = 0.0, uphi = 0.0, uz = 0.0;
    double dur_dr = 0.0, dur_dphi = 0.0, dur_dz = 0.0;
    double duphi_dr = 0.0, duphi_dphi = 0.0, duphi_dz = 0.0;
    double duz_dr = 0.0, duz_dphi = 0.0, duz_dz = 0.0;
};

// Strain order: [eps_rr, eps_phiphi, eps_zz, gamma_rphi, gamma_rz, gamma_phiz].
// Shears are engineering shear strains (twice tensorial shear).
using CylindricalStrain = std::array<double, 6>;
using ElasticityMatrix6 = std::array<std::array<double, 6>, 6>;

CylindricalStrain cylindrical_small_strain(const CylindricalDisplacementJet& u);
ElasticityMatrix6 isotropic_elasticity_3d(double young_modulus_pa, double poisson_ratio);
double elastic_strain_energy_density(const CylindricalStrain& strain,
                                     const ElasticityMatrix6& elasticity);

} // namespace BellOptimizer::Math
