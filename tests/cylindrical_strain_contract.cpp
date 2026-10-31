#include "math/cylindrical_strain.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer::Math;
static void req(bool v, const char* msg) { if (!v) throw std::runtime_error(msg); }
static void zero(const CylindricalStrain& e, const char* msg) {
    for (double x : e) req(std::abs(x) < 1e-12, msg);
}
int main() {
    // Uniform translation along z has no strain.
    CylindricalDisplacementJet uz; uz.r=0.37; uz.uz=0.4;
    zero(cylindrical_small_strain(uz), "axial rigid translation must have zero strain");
    // Rigid translation in x: ur=U cos(phi), uphi=-U sin(phi).
    const double r=0.37, U=0.02, phi=0.61;
    CylindricalDisplacementJet tx; tx.r=r;
    tx.ur=U*std::cos(phi); tx.uphi=-U*std::sin(phi);
    tx.dur_dphi=-U*std::sin(phi); tx.duphi_dphi=-U*std::cos(phi);
    zero(cylindrical_small_strain(tx), "rigid transverse translation must have zero strain");
    // Rigid rotation about z: uphi=Omega*r, d(uphi)/dr=Omega.
    CylindricalDisplacementJet rot; rot.r=r; rot.uphi=1.7*r; rot.duphi_dr=1.7;
    zero(cylindrical_small_strain(rot), "rigid axial rotation must have zero strain");
    const auto C=isotropic_elasticity_3d(110e9,0.34);
    CylindricalStrain axial{{0,0,1e-5,0,0,0}};
    const double w=elastic_strain_energy_density(axial,C);
    req(w>0.0 && std::isfinite(w), "positive strain must have positive finite energy");
    bool bad=false; try { (void)isotropic_elasticity_3d(1.0,0.5); } catch(const std::invalid_argument&) { bad=true; }
    req(bad,"incompressible singular endpoint nu=0.5 must be rejected");
    std::cout << "cylindrical_strain_contract: rigid-body invariants and isotropic energy PASS\n";
}
