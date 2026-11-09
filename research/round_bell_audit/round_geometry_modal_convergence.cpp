#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "fem/internal_modal_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr int kModes = 10;
struct Level { int n = 0; BellOptimizer::ModalData modes; };
struct PairMetric { int coarse_mode = 0, fine_mode = 0; double freq_change = 0, corr = 0, shape_error = 0; };

std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> v; std::stringstream ss(s); std::string x;
    while (std::getline(ss, x, ',')) v.push_back(x);
    return v;
}

bool buildLevel(const BellGeometry& geometry, int n, const BellMaterial& material, Level& level) {
    BellMeshConfig config; config.nNodes = n;
    BellGenerator generator(material, config); generator.setGeometry(geometry);
    if (!generator.validateGeometry().passed || !generator.buildMesh()) return false;
    const auto& mesh = generator.getMeshData();
    BellOptimizer::InternalModalSolver solver(true, kModes, true);
    BellOptimizer::Excitation::MeridionalInfluence influence;
    BellOptimizer::ModalProblem problem{mesh, material, ".", "round_geometry_modal_convergence", &influence};
    if (!solver.solve(problem, level.modes) || level.modes.modes.size() < kModes) return false;
    level.n = n; return true;
}

void sample(const BellOptimizer::ModalMode& mode, int n, double t, std::array<double,4>& out) {
    const double x = std::clamp(t, 0.0, 1.0) * static_cast<double>(n - 1);
    const int i0 = std::min(n - 1, static_cast<int>(std::floor(x)));
    const int i1 = std::min(n - 1, i0 + 1);
    const double a = x - static_cast<double>(i0);
    auto interp = [&](const std::vector<double>& v, int off) {
        const double p = v.at(static_cast<std::size_t>(off + i0));
        const double q = v.at(static_cast<std::size_t>(off + i1));
        return p + a * (q - p);
    };
    out = {interp(mode.nodal_radial,0), interp(mode.nodal_axial,0),
           interp(mode.nodal_radial,n), interp(mode.nodal_axial,n)};
}

PairMetric compare(const BellOptimizer::ModalMode& a, int na,
                   const BellOptimizer::ModalMode& b, int nb, int ia, int ib) {
    double dot=0, aa=0, bb=0; std::vector<double> va, vb;
    for (int i=0;i<nb;++i) {
        std::array<double,4> x{},y{};
        sample(a,na,static_cast<double>(i)/(nb-1),x); sample(b,nb,static_cast<double>(i)/(nb-1),y);
        va.insert(va.end(),x.begin(),x.end()); vb.insert(vb.end(),y.begin(),y.end());
    }
    for (std::size_t i=0;i<va.size();++i) { dot+=va[i]*vb[i]; aa+=va[i]*va[i]; bb+=vb[i]*vb[i]; }
    PairMetric p; p.coarse_mode=ia; p.fine_mode=ib;
    if (!(aa>0) || !(bb>0)) return p;
    p.corr=std::abs(dot)/std::sqrt(aa*bb);
    const double scale=dot/aa; double residual=0;
    for (std::size_t i=0;i<va.size();++i) { const double d=vb[i]-scale*va[i]; residual+=d*d; }
    p.shape_error=std::sqrt(residual/bb);
    const double fc=a.frequency_hz, ff=b.frequency_hz;
    p.freq_change=std::abs(ff-fc)/std::max(std::abs(ff),1e-12);
    return p;
}

int bitCount(unsigned value) {
    int count = 0;
    while (value) { count += static_cast<int>(value & 1u); value >>= 1u; }
    return count;
}

std::vector<PairMetric> assign(const Level& a, const Level& b) {
    // Exact maximum-weight one-to-one assignment for ten modes (portable C++17 bitmask DP).
    std::array<std::array<PairMetric,kModes>,kModes> metrics{};
    for(int i=0;i<kModes;++i) for(int j=0;j<kModes;++j)
        metrics[i][j]=compare(a.modes.modes[i],a.n,b.modes.modes[j],b.n,i+1,j+1);
    constexpr int states=1<<kModes;
    std::array<double,states> dp; dp.fill(-1e100); dp[0]=0;
    std::array<int,states> prev_mask{}, prev_choice{};
    for(int mask=0;mask<states;++mask) {
        const int i=__builtin_popcount(static_cast<unsigned>(mask));
        if(i>=kModes || dp[mask]<-1e90) continue;
        for(int j=0;j<kModes;++j) if(!(mask&(1<<j))) {
            const int next=mask|(1<<j); const double score=dp[mask]+metrics[i][j].corr;
            if(score>dp[next]) { dp[next]=score; prev_mask[next]=mask; prev_choice[next]=j; }
        }
    }
    std::vector<PairMetric> out(kModes); int mask=states-1;
    for(int i=kModes-1;i>=0;--i) { const int j=prev_choice[mask]; out[i]=metrics[i][j]; mask=prev_mask[mask]; }
    return out;
}
}

int main(int argc,char**argv) {
    const std::string input=argc>1?argv[1]:"round_geometry_dataset.csv";
    const std::string output=argc>2?argv[2]:"round_geometry_modal_convergence.csv";
    std::ifstream in(input); if(!in) { std::cerr<<"Cannot open dataset "<<input<<'\n'; return 1; }
    std::string line; std::getline(in,line);
    std::ofstream csv(output); if(!csv) return 2;
    csv<<"geometry_index,tag,coarse_nodes,fine_nodes,coarse_mode,fine_mode,coarse_hz,fine_hz,relative_frequency_change,shape_correlation,normalized_shape_l2_error\n";
    csv<<std::setprecision(12);
    const std::vector<int> resolutions{41,51,61,71,81};
    BellMaterial material; int geometries=0, failures=0;
    std::array<double,4> max_df{}; std::array<double,4> min_corr{1,1,1,1}; std::array<double,4> max_shape{};
    while(std::getline(in,line)) {
        if(line.empty()) continue;
        const auto f=split(line); if(f.size()<17 || f[16]!="PASS") continue;
        BellGeometry g;
        g.H=std::stod(f[2]); g.R_mouth=std::stod(f[3]); g.t_top=std::stod(f[4]); g.t_bottom=std::stod(f[5]);
        g.r_knots={0.0,std::stod(f[6]),std::stod(f[7]),std::stod(f[8]),1.0}; g.ellipticity=std::stod(f[9]);
        std::vector<Level> levels;
        bool ok=true;
        for(int n:resolutions) { Level l; if(!buildLevel(g,n,material,l)) { ok=false; break; } levels.push_back(std::move(l)); }
        if(!ok) { ++failures; std::cerr<<"MODAL_FAIL geometry="<<f[0]<<" tag="<<f[1]<<'\n'; continue; }
        ++geometries;
        for(std::size_t p=1;p<levels.size();++p) {
            const auto pairs=assign(levels[p-1],levels[p]);
            for(const auto& m:pairs) {
                const auto& cm=levels[p-1].modes.modes[m.coarse_mode-1];
                const auto& fm=levels[p].modes.modes[m.fine_mode-1];
                csv<<f[0]<<','<<f[1]<<','<<levels[p-1].n<<','<<levels[p].n<<','<<m.coarse_mode<<','<<m.fine_mode<<','
                   <<cm.frequency_hz<<','<<fm.frequency_hz<<','<<m.freq_change<<','<<m.corr<<','<<m.shape_error<<'\n';
                max_df[p-1]=std::max(max_df[p-1],m.freq_change);
                min_corr[p-1]=std::min(min_corr[p-1],m.corr);
                max_shape[p-1]=std::max(max_shape[p-1],m.shape_error);
            }
        }
        if(geometries%5==0) std::cout<<"completed_geometries="<<geometries<<'\n';
    }
    std::cout<<"modal_geometries="<<geometries<<" modal_failures="<<failures<<'\n';
    for(int p=0;p<4;++p) std::cout<<resolutions[p]<<"->"<<resolutions[p+1]
        <<" max_frequency_change_percent="<<max_df[p]*100.0
        <<" min_shape_correlation="<<min_corr[p]
        <<" max_shape_error_percent="<<max_shape[p]*100.0<<'\n';
    return geometries>=60 && failures==0 && csv ? 0 : 3;
}
