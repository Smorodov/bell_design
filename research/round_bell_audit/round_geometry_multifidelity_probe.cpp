#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/cylindrical_bell_model.hpp"
#include "fem/cylindrical_modal_evaluator.hpp"
#include "fem/internal_modal_solver.hpp"
#include "core/modal_tracking.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using namespace BellOptimizer;
namespace {
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while(std::getline(ss,x,',')) v.push_back(x);
    return v;
}
BellGeometry parseGeometry(const std::vector<std::string>& f) {
    BellGeometry g; g.H=std::stod(f[2]); g.R_mouth=std::stod(f[3]);
    g.t_top=std::stod(f[4]); g.t_bottom=std::stod(f[5]);
    g.r_knots={0.0,std::stod(f[6]),std::stod(f[7]),std::stod(f[8]),1.0};
    g.ellipticity=std::stod(f[9]); return g;
}
bool solveFast(const BellGeometry& g,const BellMaterial& material,const BellMeshConfig& mesh,ModalData& data) {
    BellGenerator generator(material,mesh); generator.setGeometry(g);
    if(!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    Excitation::MeridionalInfluence influence;
    ModalProblem problem{generator.getMeshData(),material,".","round_mf_fast",&influence};
    InternalModalSolver solver(true,8,true);
    return solver.solve(problem,data) && data.modes.size()>=5;
}
bool solveDirect(const BellGeometry& g,const BellMaterial& material,const BellMeshConfig& mesh,ModalData& data) {
    BellGenerator generator(material,mesh); generator.setGeometry(g);
    if(!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    const auto model=FEM::make_cylindrical_bell_model(g,material,mesh);
    Excitation::MeridionalInfluence provider;
    const auto influence=provider.make(model.mesh_si);
    FEM::CylindricalModalEvaluationConfig cfg;
    cfg.harmonics={0}; cfg.angular_samples=32; cfg.requested_modes=8;
    cfg.max_iterations=1600; cfg.tolerance=1e-6; cfg.minimum_frequency_hz=1.0;
    return FEM::evaluate_cylindrical_bell_modes(model,influence,cfg,data) && data.modes.size()>=5;
}
double objective(const std::array<double,5>& f) {
    if(!(f[1]>0.0)) return 1e9;
    const std::array<double,4> target{0.5,1.25,1.5,2.0};
    const std::array<int,4> index{0,2,3,4}; double score=0;
    for(std::size_t i=0;i<4;++i) {
        const double ratio=f[static_cast<std::size_t>(index[i])]/f[1];
        if(!(ratio>0.0)||!std::isfinite(ratio)) return 1e9;
        const double e=std::log(ratio/target[i]); score+=e*e;
    }
    return score;
}
std::vector<double> ranks(const std::vector<double>& x) {
    std::vector<std::size_t> order(x.size()); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return x[a]<x[b];});
    std::vector<double> r(x.size()); for(std::size_t i=0;i<order.size();++i) r[order[i]]=double(i+1); return r;
}
double spearman(const std::vector<double>& a,const std::vector<double>& b) {
    if(a.size()!=b.size()||a.size()<2) return 0;
    const auto x=ranks(a),y=ranks(b); double mx=0,my=0;
    for(double v:x)mx+=v; for(double v:y)my+=v; mx/=x.size(); my/=y.size();
    double c=0,vx=0,vy=0; for(std::size_t i=0;i<x.size();++i){double dx=x[i]-mx,dy=y[i]-my;c+=dx*dy;vx+=dx*dx;vy+=dy*dy;}
    return vx>0&&vy>0?c/std::sqrt(vx*vy):0;
}
}

int main(int argc,char**argv) {
    const std::string input=argc>1?argv[1]:"round_geometry_dataset.csv";
    const std::string output=argc>2?argv[2]:"round_geometry_multifidelity.csv";
    std::ifstream in(input); if(!in)return 1; std::string line; std::getline(in,line);
    std::ofstream csv(output); if(!csv)return 2;
    csv<<"geometry_index,tag,fast_objective,tracked_direct_objective,min_shape_similarity,mean_shape_similarity,max_first5_frequency_error,ambiguous,fast_f1_hz,direct_tracked_f1_hz\n";
    csv<<std::setprecision(12);
    BellMaterial material; BellMeshConfig mesh; mesh.nNodes=15;
    std::vector<double> fast_scores,direct_scores;
    double min_similarity=1.0,max_frequency_error=0.0,sum_mean=0.0;
    std::size_t ambiguous=0,success=0,failures=0;
    while(std::getline(in,line)) {
        if(line.empty())continue; const auto f=split(line);
        if(f.size()<17 || f[16]!="PASS")continue;
        const auto g=parseGeometry(f); ModalData fast,direct;
        if(!solveFast(g,material,mesh,fast)){++failures;std::cerr<<"FAST_FAIL "<<f[0]<<' '<<f[1]<<'\n';continue;}
        if(!solveDirect(g,material,mesh,direct)){++failures;std::cerr<<"DIRECT_FAIL "<<f[0]<<' '<<f[1]<<'\n';continue;}
        const auto tracking=trackModalFamilies(fast,direct,0.10,0.01);
        if(!tracking.used_shape_tracking||tracking.matched_modes<5){++failures;std::cerr<<"TRACK_FAIL "<<f[0]<<' '<<f[1]<<'\n';continue;}
        std::array<double,5> ff{},fd{}; double maxerr=0;
        for(std::size_t i=0;i<5;++i){ff[i]=fast.modes[i].frequency_hz;fd[i]=direct.modes[i].frequency_hz;
            maxerr=std::max(maxerr,std::abs(ff[i]-fd[i])/std::max(std::abs(fd[i]),1e-12));}
        const double fs=objective(ff),ds=objective(fd); const bool amb=tracking.ambiguous||tracking.minimum_similarity<0.75;
        fast_scores.push_back(fs);direct_scores.push_back(ds);min_similarity=std::min(min_similarity,tracking.minimum_similarity);
        max_frequency_error=std::max(max_frequency_error,maxerr);sum_mean+=tracking.mean_similarity;
        if(amb)++ambiguous;++success;
        csv<<f[0]<<','<<f[1]<<','<<fs<<','<<ds<<','<<tracking.minimum_similarity<<','<<tracking.mean_similarity
           <<','<<maxerr<<','<<(amb?1:0)<<','<<ff[0]<<','<<fd[0]<<'\n';
        if(success%10==0)std::cout<<"completed="<<success<<'\n';
    }
    const double rho=spearman(fast_scores,direct_scores);
    std::cout<<"accepted_geometries="<<success<<" failures="<<failures<<" ambiguous="<<ambiguous<<'\n'
             <<"spearman_objective_rank="<<rho<<" max_first5_frequency_error="<<max_frequency_error
             <<" min_shape_similarity="<<min_similarity<<" mean_candidate_mean_similarity="
             <<(success?sum_mean/success:0)<<'\n'
             <<"scope=m0_only; same_15_node_profile_mesh; boundary_trace_modal_tracking; no_full_volume_mass_metric\n";
    return success>=60 && failures==0 && std::isfinite(rho) && csv?0:3;
}
