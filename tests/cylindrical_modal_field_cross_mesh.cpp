#include "fem/cylindrical_modal_field.hpp"
#include "fem/cylindrical_global_assembly.hpp"
#include "fem/cylindrical_global_modes.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

static void require(bool value,const char* message) { if(!value) throw std::runtime_error(message); }

static BellMeshData make_patch(int nr,int nz,double zmax=0.10) {
    BellMeshData mesh;
    for(int j=0;j<=nz;++j) for(int i=0;i<=nr;++i) {
        const double r=0.20+0.10*static_cast<double>(i)/nr;
        const double z=zmax*static_cast<double>(j)/nz;
        mesh.nodes.push_back({static_cast<int>(mesh.nodes.size())+1,r,z});
    }
    const int stride=nr+1;
    for(int j=0;j<nz;++j) for(int i=0;i<nr;++i) {
        const int n1=j*stride+i+1,n2=n1+1,n4=(j+1)*stride+i+1,n3=n4+1;
        mesh.elements.push_back({n1,n2,n3,n4});
    }
    return mesh;
}

static std::vector<int> clamp_lower_edge(const BellMeshData& mesh,std::size_t dpn) {
    std::vector<int> fixed;
    for(std::size_t i=0;i<mesh.nodes.size();++i) if(std::abs(mesh.nodes[i].z)<1e-12)
        for(std::size_t d=0;d<dpn;++d) fixed.push_back(static_cast<int>(i*dpn+d));
    return fixed;
}

static std::vector<double> make_mode(const BellMeshData& mesh,std::size_t basis_count,
                                     std::size_t component,double scale) {
    std::vector<double> v(mesh.nodes.size()*3u*basis_count,0.0);
    for(std::size_t i=0;i<mesh.nodes.size();++i) {
        const auto& p=mesh.nodes[i];
        v[(i*3u+component)*basis_count]=scale*(1.0+2.0*p.r+3.0*p.z);
    }
    return v;
}

int main() {
    const auto coarse=make_patch(2,2);
    const auto fine=make_patch(5,4);
    const auto clipped=make_patch(5,4,0.080);
    const std::vector<int> h0{0}, h02{0,2};
    auto a0=make_mode(coarse,1,0,1.0);
    auto b0=make_mode(fine,3,0,1.0);
    auto a1=make_mode(coarse,1,2,0.7);
    auto b1=make_mode(fine,3,2,0.7);
    const CylindricalModalField fa0{&coarse,&h0,&a0,0.85,7800.0};
    const CylindricalModalField fb0{&fine,&h02,&b0,0.85,7800.0};
    const auto overlap=compare_cylindrical_modal_fields(fa0,fb0,64,0.98);
    require(overlap.valid,"same physical field on different meshes/bases was rejected");
    require(overlap.normalized_overlap>0.999999,"cross-mesh mass overlap is not close to unity");
    require(overlap.common_volume_fraction>0.999999,"identical domains do not have full common-volume coverage");

    const CylindricalModalField fa1{&coarse,&h0,&a1,0.85,7800.0};
    const CylindricalModalField fb1{&fine,&h02,&b1,0.85,7800.0};
    double coverage=0.0;
    const double pair=cylindrical_modal_subspace_similarity({fa0,fa1},{fb0,fb1},64,0.98,&coverage);
    require(pair>0.999999,"cross-mesh two-dimensional modal subspace similarity is not unity");
    require(coverage>0.999999,"cross-mesh subspace coverage is incomplete");

    auto clipped_mode=make_mode(clipped,3,0,1.0);
    const CylindricalModalField fclip{&clipped,&h02,&clipped_mode,0.85,7800.0};
    const auto mismatch=compare_cylindrical_modal_fields(fa0,fclip,64,0.98);
    require(!mismatch.valid,"different physical domains passed the common-volume gate");
    require(mismatch.common_volume_fraction<0.98,"domain mismatch was not detected by coverage");

    bool alias_rejected=false;
    const std::vector<int> high_harmonics{0,8};
    auto high_mode=make_mode(fine,3,0,1.0);
    const CylindricalModalField high_field{&fine,&high_harmonics,&high_mode,0.85,7800.0};
    try { (void)compare_cylindrical_modal_fields(high_field,high_field,16,0.98); }
    catch(const std::invalid_argument&) { alias_rejected=true; }
    require(alias_rejected,"angular quadrature aliasing was not rejected");
    bool geometry_mismatch_rejected=false;
    const CylindricalModalField wrong_q{&fine,&h02,&b0,0.90,7800.0};
    try { (void)compare_cylindrical_modal_fields(fa0,wrong_q,64,0.98); }
    catch(const std::invalid_argument&) { geometry_mismatch_rejected=true; }
    require(geometry_mismatch_rejected,"different ellipticity values were accepted in one physical mass metric");

    // A real FEM comparison on two nested meshes. The mode number is not
    // assumed to be identical: choose the best physical-field continuation
    // among the first four fine-mesh candidates.
    const auto fem_coarse=make_patch(4,4);
    const auto fem_fine=make_patch(8,8);
    const std::vector<int> fem_harmonics{0};
    constexpr double young=1.0e11, poisson=0.30, density=7800.0, q=0.85;
    const auto coarse_system=assemble_cylindrical_global_system(fem_coarse,young,poisson,density,fem_harmonics,64,0.0,q);
    const auto fine_system=assemble_cylindrical_global_system(fem_fine,young,poisson,density,fem_harmonics,64,0.0,q);
    const auto coarse_modes=solve_cylindrical_global_modes(coarse_system,
        clamp_lower_edge(fem_coarse,coarse_system.dofs_per_node),4,1200,1e-4,true);
    const auto fine_modes=solve_cylindrical_global_modes(fine_system,
        clamp_lower_edge(fem_fine,fine_system.dofs_per_node),4,1200,1e-4,true);
    require(coarse_modes.converged && fine_modes.converged,"cross-mesh FEM eigensolve failed");
    double best_fem_overlap=0.0; std::size_t best_coarse=0,best_fine=0;
    for(std::size_t i=0;i<coarse_modes.mode_shapes.size();++i)
        for(std::size_t j=0;j<fine_modes.mode_shapes.size();++j) {
            const CylindricalModalField ca{&fem_coarse,&fem_harmonics,&coarse_modes.mode_shapes[i],q,density};
            const CylindricalModalField fb{&fem_fine,&fem_harmonics,&fine_modes.mode_shapes[j],q,density};
            const auto ov=compare_cylindrical_modal_fields(ca,fb,64,0.98);
            if(ov.valid && ov.normalized_overlap>best_fem_overlap) {
                best_fem_overlap=ov.normalized_overlap; best_coarse=i; best_fine=j;
            }
        }
    require(best_fem_overlap>0.70,"real FEM modal field did not continue across mesh refinement");

    std::cout<<"cross_mesh_overlap="<<overlap.normalized_overlap
             <<" common_volume="<<overlap.common_volume_fraction
             <<" subspace="<<pair<<" clipped_coverage="<<mismatch.common_volume_fraction
             <<" fem_best_overlap="<<best_fem_overlap<<" fem_modes="<<best_coarse<<","<<best_fine<<'\n';
    std::cout<<"cylindrical_modal_field_cross_mesh: PASS\n";
    return 0;
}
