#include "fem/cylindrical_q4_element.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main() {
    using namespace BellOptimizer::FEM;
    const std::array<MeridionalNode,4> q{{{0.20,0.00},{0.30,0.00},{0.30,0.10},{0.20,0.10}}};
    const auto e=assemble_cylindrical_q4(q, 1.0e11, 0.30, 7800.0, {0,1,2}, 128);
    const std::size_t n=e.mass.size();
    require(n==12u*5u,"expected three cylindrical components times five real Fourier basis functions per node");
    double kscale=0.0, masym=0.0, kasym=0.0, mdiag=1e300;
    for(std::size_t i=0;i<n;++i) {
        require(e.mass[i][i]>0.0,"mass diagonal must be positive");
        mdiag=std::min(mdiag,e.mass[i][i]);
        for(std::size_t j=0;j<n;++j) {
            kscale=std::max(kscale,std::abs(e.stiffness[i][j]));
            kasym=std::max(kasym,std::abs(e.stiffness[i][j]-e.stiffness[j][i]));
            masym=std::max(masym,std::abs(e.mass[i][j]-e.mass[j][i]));
        }
    }
    require(kasym/std::max(kscale,1e-300)<1e-12,"stiffness must be symmetric");
    require(masym<1e-12,"mass must be symmetric");

    // Rigid translation in global x: u_r=cos(phi), u_phi=-sin(phi), u_z=0.
    // For the {0,1,2} basis, cos(1 phi) is basis index 1 and sin(1 phi) index 2.
    std::vector<double> x(n,0.0);
    const std::size_t nb=5;
    for(std::size_t node=0;node<4;++node) {
        x[(node*3u+0u)*nb+1u]=1.0;
        x[(node*3u+1u)*nb+2u]=-1.0;
    }
    double energy=0.0;
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j) energy+=x[i]*e.stiffness[i][j]*x[j];
    require(std::abs(energy)/std::max(kscale,1.0)<1e-11,"global rigid translation must have near-zero elastic energy");

    const auto elliptical = assemble_cylindrical_q4(q, 1.0e11, 0.30, 7800.0, {0,1,2}, 128, 0.85);
    double elliptical_scale=0.0, elliptical_energy=0.0;
    for(std::size_t i=0;i<n;++i) {
        elliptical_scale=std::max(elliptical_scale,std::abs(elliptical.stiffness[i][i]));
        for(std::size_t j=0;j<n;++j) elliptical_energy += x[i]*elliptical.stiffness[i][j]*x[j];
    }
    require(std::abs(elliptical_energy)/std::max(elliptical_scale,1.0)<1e-11,
            "global rigid translation must remain strain-free under elliptical geometry mapping");

    bool bad=false;
    try { (void)assemble_cylindrical_q4(q,1e11,0.3,0.0,{0,1,2}); }
    catch(const std::invalid_argument&) { bad=true; }
    require(bad,"zero density must be rejected");
    bad=false;
    try { (void)assemble_cylindrical_q4(q,1e11,0.3,7800.0,{0,1,2},64,1.2); }
    catch(const std::invalid_argument&) { bad=true; }
    require(bad,"ellipticity above one must be rejected");
    std::cout << "cylindrical_q4_element_contract: symmetric K/M, circular/elliptical rigid translation and input validation PASS\n";
}
