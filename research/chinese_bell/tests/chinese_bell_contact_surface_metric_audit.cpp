#include "fem/cylindrical_bell_model.hpp"
#include "excitation/chinese_bell_strike.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;
namespace {
constexpr double PI = 3.1415926535897932384626433832795;
void req(bool x,const char* m){if(!x)throw std::runtime_error(m);}
BellGeometry geom(double q){BellGeometry g;g.H=658.3175709034034;g.R_mouth=355.9758748202059;g.t_top=22.85675347179014;g.t_bottom=46.62524667849419;g.r_knots={0.0,0.130778,0.216866,0.45905,1.0};g.ellipticity=q;return g;}
struct Metric {double meridional=0.0;double azimuthal=0.0;};
Metric exact_metric(const std::vector<Point2D>& p,std::size_t i,double q,double phi){
    const std::size_t j0=i==0?0:i-1, j1=std::min(i+1,p.size()-1);
    const double dr=p[j1].r-p[j0].r, dz=p[j1].z-p[j0].z;
    const double drdz=std::abs(dz)>1e-14?dr/dz:0.0;
    const double c=std::cos(phi),s=std::sin(phi), invq2=1.0/(q*q);
    const double d=c*c+s*s*invq2, alpha=1.0/std::sqrt(d);
    const double dlog=-s*c*(invq2-1.0)/d;
    return {std::sqrt(1.0+alpha*alpha*drdz*drdz), alpha*std::max(0.0,p[i].r)*std::sqrt(1.0+dlog*dlog)};
}
std::vector<double> load(const BellMeshData&m,const std::vector<int>&h,double q,double z,double az,double mw,double aw){
    Excitation::ChineseBellStrikeConfig c;c.center_z=z;c.azimuth_rad=az;c.meridional_width=mw;c.angular_width_rad=aw;
    return Excitation::make_chinese_bell_strike_influence(m,h,q,c,128).dof_values;
}
double rel(const std::vector<double>&a,const std::vector<double>&b){req(a.size()==b.size(),"size mismatch");double n=0,d=0;for(std::size_t i=0;i<a.size();++i){const double x=a[i]-b[i];n+=x*x;d+=b[i]*b[i];}return std::sqrt(n/std::max(d,1e-300));}
}
int main(){try{
    BellMaterial mat; BellMeshConfig mc; mc.nNodes=35; const std::vector<int> h={0,2,4}; const double q=.70;
    auto model=make_cylindrical_bell_model(geom(q),mat,mc);
    const double f=.70; const double z=model.mesh_si.outer_profile.front().z+f*(model.mesh_si.outer_profile.back().z-model.mesh_si.outer_profile.front().z);
    auto it=std::min_element(model.mesh_si.outer_profile.begin(),model.mesh_si.outer_profile.end(),[&](const auto&a,const auto&b){return std::abs(a.z-z)<std::abs(b.z-z);});
    const std::size_t i=static_cast<std::size_t>(it-model.mesh_si.outer_profile.begin());
    std::cout<<std::setprecision(12)<<"width_m,azimuth,exact_meridional_metric,exact_azimuthal_metric,approx_azimuthal_metric,legacy_map_error,exact_map_error\n";
    double max_legacy=0.0,max_exact=0.0;
    for(double w:{.020,.035,.050}) for(double az:{0.0,0.5*PI}){
        const auto em=exact_metric(model.mesh_si.outer_profile,i,q,az);
        const double approx_m=std::sqrt(1.0 + std::pow((model.mesh_si.outer_profile[std::min(i+1,model.mesh_si.outer_profile.size()-1)].r-model.mesh_si.outer_profile[i==0?0:i-1].r)/std::max(1e-14,model.mesh_si.outer_profile[std::min(i+1,model.mesh_si.outer_profile.size()-1)].z-model.mesh_si.outer_profile[i==0?0:i-1].z),2));
        const double approx_a=std::max(0.0,it->r);
        const double exact_mw=w/em.meridional, exact_aw=w/em.azimuthal;
        const double legacy_mw=w/approx_m, legacy_aw=w/approx_a;
        auto ref=load(model.mesh_si,h,q,z,az,exact_mw,exact_aw);
        auto legacy=load(model.mesh_si,h,q,z,az,legacy_mw,legacy_aw);
        auto exact=load(model.mesh_si,h,q,z,az,exact_mw,exact_aw);
        const double el=rel(legacy,ref), ee=rel(exact,ref);
        max_legacy=std::max(max_legacy,el); max_exact=std::max(max_exact,ee);
        std::cout<<w<<','<<az<<','<<em.meridional<<','<<em.azimuthal<<','<<approx_a<<','<<el<<','<<ee<<'\n';
        req(std::isfinite(el)&&std::isfinite(ee),"non-finite metric audit");
        req(ee < 1e-12,"exact physical mapping is not self-consistent");
    }
    req(max_legacy>1e-4,"legacy approximation unexpectedly exact");
    std::cout<<"chinese_bell_contact_surface_metric_audit: PASS\n";
    std::cout<<"max_legacy_mapping_error="<<max_legacy<<"\n";
    std::cout<<"exact_mapping_error="<<max_exact<<"\n";
    std::cout<<"status=diagnostic_surface_metric_correction_no_physical_threshold\n";
    return 0;
}catch(const std::exception&e){std::cerr<<"chinese_bell_contact_surface_metric_audit: FAIL: "<<e.what()<<'\n';return 1;}}
