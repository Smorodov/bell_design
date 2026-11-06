#include "acoustics/acoustic_field_cache.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::Acoustics;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
ModalMode makeMode(double frequency,double scale){ModalMode m;m.frequency_hz=frequency;m.outer_surface_harmonics={0,1,-1,2,-2};const std::size_t points=5,basis=m.outer_surface_harmonics.size();m.outer_surface_fourier_radial.resize(points*basis);m.outer_surface_fourier_axial.resize(points*basis);for(std::size_t p=0;p<points;++p)for(std::size_t b=0;b<basis;++b){const double x=0.17*static_cast<double>((p+1)*(b+2));m.outer_surface_fourier_radial[p*basis+b]=scale*0.0006*std::sin(x);m.outer_surface_fourier_axial[p*basis+b]=scale*0.0003*std::cos(x*0.8);}return m;}
}
int main(){const std::string path="acoustic_field_cache_contract.bin";try{
 const auto m1=makeMode(640,1),m2=makeMode(913,0.7);const std::vector<const ModalMode*> modes{&m1,&m2};const std::vector<Point2D> profile{{35,0},{60,90},{95,190},{135,300},{180,420}};const AcousticFieldCacheIdentity id{"geometry-sha256-demo","bronze-preset-v2","fem-modal-v4"};
 const auto built=buildAcousticFieldCache(modes,profile,id,91);require(built.valid(),"built cache invalid");saveAcousticFieldCache(path,built);const auto loaded=loadAcousticFieldCache(path,id);require(loaded.valid()&&loaded.modes.size()==2,"loaded cache invalid or mode count changed");
 // The same geometry-dependent field must serve several receiver distances,
 // centre polar/azimuth positions, head yaws and both ears without rebuilding.
 for(const auto& receiver:std::vector<std::array<double,4>>{{3.0,0.82,-0.31,0.18},{1.0,1.20,0.74,0.18},{6.0,0.45,-1.10,0.18}}){
  const auto prepared=prepareStereoRadiationFromCache(loaded,receiver[0],receiver[1],receiver[2],receiver[3]);
  require(prepared.modeCount()==loaded.modes.size(),"prepared cache mode count mismatch");
  for(std::size_t m=0;m<prepared.modeCount();++m)require(std::abs(prepared.modeFrequencyHz(m)-loaded.modes[m].frequency_hz)<1e-12,"cached frequency changed");
  for(double yaw:{-1.5707963267948966,-0.7,0.0,0.8,1.5707963267948966}){const auto transfers=prepared.query(yaw,0.23);require(transfers.size()==2,"cached transfer count mismatch");for(const auto& t:transfers)require(std::isfinite(t.left.real())&&std::isfinite(t.right.imag())&&std::isfinite(t.left_delay_s)&&std::isfinite(t.right_delay_s),"non-finite cached receiver result");}
  require(prepared.surfaceIntegralBuildCount()==0,"loading/querying cache unexpectedly rebuilt surface integrals");
 }
 bool mismatch=false;try{(void)loadAcousticFieldCache(path,{"other-geometry",id.material_fingerprint,id.solver_fingerprint});}catch(const std::runtime_error&){mismatch=true;}require(mismatch,"cache accepted a different geometry identity");
 {std::fstream f(path,std::ios::binary|std::ios::in|std::ios::out);f.seekp(-1,std::ios::end);char c=0;f.write(&c,1);}bool corrupt=false;try{(void)loadAcousticFieldCache(path,id);}catch(const std::runtime_error& e){corrupt=std::string(e.what()).find("checksum")!=std::string::npos;}require(corrupt,"cache corruption not detected");std::remove(path.c_str());
 std::cout<<"acoustic_field_cache_contract: PASS modes="<<loaded.modes.size()<<" polar_nodes_per_mode="<<loaded.modes.front().directional_field.fields.size()<<" reusable_receiver_configurations=3\n";return 0;
 }catch(const std::exception& e){std::remove(path.c_str());std::cerr<<"acoustic_field_cache_contract: FAIL: "<<e.what()<<"\n";return 1;}}
