#include "fem/cylindrical_global_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

static void require(bool v, const char* m) { if (!v) throw std::runtime_error(m); }
int main() {
    using namespace BellOptimizer;
    BellMeshData mesh;
    mesh.nodes = {{1,0.20,0.00},{2,0.30,0.00},{3,0.30,0.10},{4,0.20,0.10}};
    mesh.elements = {{1,2,3,4}};
    const auto a = FEM::assemble_cylindrical_global_system(mesh, 1.0e11, 0.30, 7800.0, {0,1,2}, 32);
    const auto b = FEM::assemble_cylindrical_global_system(mesh, 1.0e11, 0.30, 7800.0, {0,1,2}, 64);
    require(a.stiffness.n == 60 && a.mass.n == 60, "global DOF dimension incorrect");
    require(a.stiffness.values.size() < static_cast<std::size_t>(a.stiffness.n)*a.stiffness.n,
            "global stiffness should be stored sparsely");
    double kscale=0.0, asym=0.0, mass_asym=0.0, conv=0.0, mscale=0.0;
    for (int i=0;i<a.stiffness.n;++i) for (int j=0;j<a.stiffness.n;++j) {
        const double kij=a.stiffness.get(i,j), mij=a.mass.get(i,j);
        kscale=std::max(kscale,std::abs(kij)); mscale=std::max(mscale,std::abs(mij));
        asym=std::max(asym,std::abs(kij-a.stiffness.get(j,i)));
        mass_asym=std::max(mass_asym,std::abs(mij-a.mass.get(j,i)));
        conv=std::max(conv,std::abs(kij-b.stiffness.get(i,j)));
        require(mij >= -1e-12*std::max(mscale,1.0), "unexpected negative mass coefficient");
    }
    require(asym/std::max(kscale,1.0)<1e-12, "global stiffness is not symmetric");
    require(mass_asym/std::max(mscale,1.0)<1e-12, "global mass is not symmetric");
    require(conv/std::max(kscale,1.0)<1e-11, "angular quadrature is not converged for the selected Fourier basis");
    // Uniform global x translation: ur=cos(phi), uphi=-sin(phi) at every node.
    std::vector<double> rigid(static_cast<std::size_t>(a.stiffness.n),0.0);
    const std::size_t nb=a.basis_count;
    for (std::size_t node=0;node<mesh.nodes.size();++node) {
        rigid[node*a.dofs_per_node + 0u*nb + 1u] = 1.0;
        rigid[node*a.dofs_per_node + 1u*nb + 2u] = -1.0;
    }
    const auto kr=a.stiffness.multiply(rigid);
    double energy=0.0; for (std::size_t i=0;i<rigid.size();++i) energy+=rigid[i]*kr[i];
    require(std::abs(energy)/std::max(kscale,1.0)<1e-11, "global rigid translation has nonzero energy");
    bool bad=false; mesh.elements[0].n4=999;
    try { (void)FEM::assemble_cylindrical_global_system(mesh,1e11,0.3,7800.0,{0,1},32); }
    catch(const std::out_of_range&) { bad=true; }
    require(bad,"invalid element connectivity was not rejected");
    std::cout << "cylindrical_global_assembly_contract: sparse assembly, symmetry, rigid motion, quadrature convergence PASS\n";
}
