#include "bell_generator.h"
#include "fem/internal_modal_solver.hpp"
#include "excitation/localized_strike.hpp"
#include "excitation/meridional_influence.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
struct Spectrum { double rms_log_ratio = 0.0; double cosine = 0.0; };

Spectrum compare(const std::vector<double>& a0, const std::vector<double>& b0) {
    assert(a0.size() == b0.size() && !a0.empty());
    constexpr double eps = 1e-14;
    double na=0.0, nb=0.0;
    for (size_t i=0;i<a0.size();++i) { na=std::max(na,std::abs(a0[i])); nb=std::max(nb,std::abs(b0[i])); }
    assert(na>eps && nb>eps);
    double sq=0.0,dot=0.0,aa=0.0,bb=0.0;
    for (size_t i=0;i<a0.size();++i) {
        const double a=std::abs(a0[i])/na, b=std::abs(b0[i])/nb;
        const double lr=std::log(std::max(a,eps)/std::max(b,eps)); sq+=lr*lr;
        dot+=a*b; aa+=a*a; bb+=b*b;
    }
    return {std::sqrt(sq/a0.size()), dot/std::sqrt(aa*bb)};
}
}

int main() {
    BellMaterial material{1.1e4,0.34,8.6e-9};
    BellGeometry geometry;
    geometry.H=658.3175709034034; geometry.R_mouth=355.9758748202059;
    geometry.t_top=22.85675347179014; geometry.t_bottom=46.62524667849419;
    geometry.r_knots={0.0,0.130778,0.216866,0.45905,1.0};
    BellMeshConfig cfg; cfg.nNodes=15;
    BellGenerator generator(material,cfg); generator.setGeometry(geometry);
    GeometryValidationResult validation;
    if(!generator.buildMesh(&validation)){std::cerr<<validation.reason<<"\n";return 2;}
    BellOptimizer::InternalModalSolver solver(true,9,true);
    BellOptimizer::ModalData data;
    BellOptimizer::Excitation::MeridionalInfluence influence;
    if(!solver.solve({generator.getMeshData(),material,".","localized_strike_spectrum",&influence},data) || data.modes.size()<5) return 3;

    std::vector<double> rz;
    for(const auto& m:data.modes) rz.push_back(m.part_rz*m.acoustic_volume_displacement*(2.0*3.14159265358979323846*m.frequency_hz));

    const std::vector<size_t> centers={2,5,8,11,13};
    const std::vector<double> widths={5.0,10.0,20.0,40.0};
    double best=1e300; double best_cos=0.0; size_t best_i=0; double best_w=0.0;
    for(const size_t i:centers) {
        const double z=generator.getOuterProfile()[i].z;
        for(const double w:widths) {
            BellOptimizer::Excitation::LocalizedStrikeConfig c;
            c.center_z=z; c.width=w; c.normal_fraction=1.0; c.tangent_fraction=0.0;
            const auto imp=BellOptimizer::Excitation::modalImpulse(generator.getMeshData(),data.modes,c);
            std::vector<double> pressure;
            for(size_t m=0;m<imp.size();++m) pressure.push_back(imp[m]*data.modes[m].acoustic_volume_displacement*(2.0*3.14159265358979323846*data.modes[m].frequency_hz));
            const auto s=compare(rz,pressure);
            std::cout<<"z_mm="<<std::setprecision(6)<<z<<" width_mm="<<w
                     <<" cosine="<<s.cosine<<" rms_log_ratio="<<s.rms_log_ratio<<"\n";
            if(s.rms_log_ratio<best){best=s.rms_log_ratio;best_cos=s.cosine;best_i=i;best_w=w;}
        }
    }
    assert(std::isfinite(best));
    assert(best>1e-8);
    std::cout<<"best_finite_patch_z_mm="<<generator.getOuterProfile()[best_i].z
             <<" width_mm="<<best_w<<" cosine="<<best_cos
             <<" rms_log_ratio="<<best<<"\n";
    std::cout<<"finite localized strike audit: PASS\n";
    return 0;
}
