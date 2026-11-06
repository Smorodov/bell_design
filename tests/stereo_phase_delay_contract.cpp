#include "audio/modal_pressure_synthesis.hpp"
#include "acoustics/listener_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
}
int main(){try{
 constexpr double pi=3.14159265358979323846,c=343.0,fs=48000.0;
 const auto ears=Acoustics::makeStereoListenerGeometry(2.0,0.9,0.2,0.18,0.35,c);
 const double expectedLeft=(ears.left.distance_m-2.0)/c;
 const double expectedRight=(ears.right.distance_m-2.0)/c;
 require(std::abs(ears.left.propagation_delay_s-expectedLeft)<1e-15,"left propagation delay does not match geometric path");
 require(std::abs(ears.right.propagation_delay_s-expectedRight)<1e-15,"right propagation delay does not match geometric path");
 require(std::abs((ears.right.propagation_delay_s-ears.left.propagation_delay_s)-(ears.right.distance_m-ears.left.distance_m)/c)<1e-15,"interaural time difference is inconsistent");
 Audio::ModalPressureTerm term;term.omega_rad_s=2*pi*731.0;term.alpha_per_second=0.8;term.displacement_at_release=0.003;term.velocity_at_release=-0.21;term.radiation_transfer={0.73,-0.42};const std::vector<Audio::ModalPressureTerm> modes{term};
 const double release=0.0573;const double tBefore=release+std::min(expectedLeft,expectedRight)-1.0/fs;
 require(Audio::evaluateModalPressureAt(modes,tBefore,release,std::min(expectedLeft,expectedRight))==0.0,"signal appeared before propagation arrival");
 const double local=0.00437;const double atLeft=Audio::evaluateModalPressureAt(modes,release+expectedLeft+local,release,expectedLeft);
 const double atRight=Audio::evaluateModalPressureAt(modes,release+expectedRight+local,release,expectedRight);
 const double reference=Audio::evaluateModalPressureAt(modes,release+local,release,0.0);
 require(std::abs(atLeft-reference)<1e-13*std::max(1.0,std::abs(reference)),"left channel is not a pure once-only time shift");
 require(std::abs(atRight-reference)<1e-13*std::max(1.0,std::abs(reference)),"right channel is not a pure once-only time shift");
 // The first nonzero sample is quantized to the first sample at/after arrival.
 const auto firstSample=[=](double arrival){return static_cast<long long>(std::ceil(arrival*fs-1e-12));};
 const auto leftSample=firstSample(release+expectedLeft),rightSample=firstSample(release+expectedRight);
 require(std::abs((leftSample-rightSample)/fs-(expectedLeft-expectedRight))<=1.0/fs+1e-12,"sample-grid ITD exceeds one sample quantization");
 bool rejected=false;try{(void)Audio::evaluateModalPressureAt(modes,std::numeric_limits<double>::quiet_NaN(),0.0,0.0);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"non-finite sample time was accepted");
 std::cout<<"stereo_phase_delay_contract: PASS left_delay="<<expectedLeft<<" right_delay="<<expectedRight<<" ITD_samples="<<leftSample-rightSample<<"\n";return 0;
 }catch(const std::exception& e){std::cerr<<"stereo_phase_delay_contract: FAIL: "<<e.what()<<"\n";return 1;}}
