#include "acoustic_field_cache.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace BellOptimizer::Acoustics {
namespace {
constexpr char kMagic[8] = {'B','E','L','L','A','C','F','1'};
constexpr std::uint32_t kFormatVersion = 1;
constexpr std::uint64_t kMaxPayload = 512ULL * 1024ULL * 1024ULL;
constexpr std::uint32_t kMaxString = 1024U * 1024U;
constexpr std::uint32_t kMaxModes = 10000U;
constexpr std::uint32_t kMaxNodes = 721U;
constexpr std::uint32_t kMaxCoefficients = 4097U;
constexpr double kPi = 3.14159265358979323846;

void putU32(std::vector<std::uint8_t>& b,std::uint32_t v){for(unsigned i=0;i<4;++i)b.push_back(static_cast<std::uint8_t>(v>>(8*i)));}
void putU64(std::vector<std::uint8_t>& b,std::uint64_t v){for(unsigned i=0;i<8;++i)b.push_back(static_cast<std::uint8_t>(v>>(8*i)));}
void putDouble(std::vector<std::uint8_t>& b,double v){std::uint64_t bits=0;std::memcpy(&bits,&v,8);putU64(b,bits);}
void putString(std::vector<std::uint8_t>& b,const std::string& s){if(s.empty()||s.size()>kMaxString)throw std::invalid_argument("cache fingerprints must be non-empty and bounded");putU32(b,static_cast<std::uint32_t>(s.size()));b.insert(b.end(),s.begin(),s.end());}
struct Reader{
 const std::vector<std::uint8_t>& b;std::size_t p=0;
 void need(std::size_t n){if(n>b.size()-p)throw std::runtime_error("truncated acoustic field cache");}
 std::uint32_t u32(){need(4);std::uint32_t v=0;for(unsigned i=0;i<4;++i)v|=std::uint32_t(b[p++])<<(8*i);return v;}
 std::uint64_t u64(){need(8);std::uint64_t v=0;for(unsigned i=0;i<8;++i)v|=std::uint64_t(b[p++])<<(8*i);return v;}
 double real(){auto bits=u64();double v=0;std::memcpy(&v,&bits,8);if(!std::isfinite(v))throw std::runtime_error("non-finite value in acoustic field cache");return v;}
 std::string str(){auto n=u32();if(n==0||n>kMaxString)throw std::runtime_error("invalid cache fingerprint length");need(n);std::string s(b.begin()+p,b.begin()+p+n);p+=n;return s;}
};
std::uint64_t checksum(const std::vector<std::uint8_t>& b){std::uint64_t h=14695981039346656037ULL;for(auto c:b){h^=c;h*=1099511628211ULL;}return h;}
void validateIdentity(const AcousticFieldCacheIdentity& i){if(i.geometry_fingerprint.empty()||i.material_fingerprint.empty()||i.solver_fingerprint.empty())throw std::invalid_argument("all acoustic cache fingerprints are required");}
}

bool AcousticFieldCache::valid() const {
 if(identity.geometry_fingerprint.empty()||identity.material_fingerprint.empty()||identity.solver_fingerprint.empty()||modes.empty())return false;
 for(const auto& m:modes)if(!(std::isfinite(m.frequency_hz)&&m.frequency_hz>0.0)||!m.directional_field.valid()||m.directional_field.polar_angles_rad.front()>1e-12||m.directional_field.polar_angles_rad.back()<kPi-1e-12)return false;
 return true;
}

AcousticFieldCache buildAcousticFieldCache(const std::vector<const ModalMode*>& modes,
 const std::vector<Point2D>& profile,const AcousticFieldCacheIdentity& identity,unsigned polarSamples){
 validateIdentity(identity);if(modes.empty()||profile.size()<2||polarSamples<3||polarSamples>kMaxNodes)throw std::invalid_argument("invalid acoustic field cache build input");
 AcousticFieldCache cache;cache.identity=identity;cache.modes.reserve(modes.size());
 for(const auto* mode:modes){if(!mode)throw std::invalid_argument("null mode in acoustic field cache build");CachedModalRadiationField entry;entry.frequency_hz=mode->frequency_hz;entry.directional_field=buildDirectionalRadiationPolarField(*mode,profile,1.0,polarSamples);cache.modes.push_back(std::move(entry));}
 if(!cache.valid()) throw std::runtime_error("constructed invalid acoustic field cache");
 return cache;
}

void saveAcousticFieldCache(const std::string& path,const AcousticFieldCache& cache){
 if(path.empty()||!cache.valid())throw std::invalid_argument("invalid path or acoustic field cache");
 std::vector<std::uint8_t> payload;putU32(payload,kFormatVersion);putString(payload,cache.identity.geometry_fingerprint);putString(payload,cache.identity.material_fingerprint);putString(payload,cache.identity.solver_fingerprint);putU32(payload,static_cast<std::uint32_t>(cache.modes.size()));
 for(const auto& mode:cache.modes){const auto& pf=mode.directional_field;putDouble(payload,mode.frequency_hz);putU32(payload,static_cast<std::uint32_t>(pf.fields.size()));
  for(std::size_t n=0;n<pf.fields.size();++n){const auto& f=pf.fields[n];if(f.orders.empty()||f.orders.size()!=f.coefficients.size()||f.orders.size()>kMaxCoefficients)throw std::invalid_argument("invalid harmonic coefficients in acoustic cache");putDouble(payload,pf.polar_angles_rad[n]);putDouble(payload,f.observation_radius_m);putDouble(payload,f.observation_polar_rad);putU32(payload,static_cast<std::uint32_t>(f.orders.size()));
   for(std::size_t k=0;k<f.orders.size();++k){putU32(payload,static_cast<std::uint32_t>(f.orders[k]));putDouble(payload,f.coefficients[k].real());putDouble(payload,f.coefficients[k].imag());}}
 }
 if(payload.size()>kMaxPayload) throw std::runtime_error("acoustic field cache exceeds size limit");
 std::ofstream out(path,std::ios::binary|std::ios::trunc);if(!out)throw std::runtime_error("cannot open acoustic field cache for writing");out.write(kMagic,8);std::vector<std::uint8_t> header;putU64(header,payload.size());putU64(header,checksum(payload));out.write(reinterpret_cast<const char*>(header.data()),16);out.write(reinterpret_cast<const char*>(payload.data()),static_cast<std::streamsize>(payload.size()));if(!out)throw std::runtime_error("failed writing acoustic field cache");
}

AcousticFieldCache loadAcousticFieldCache(const std::string& path,const AcousticFieldCacheIdentity& expected){
 validateIdentity(expected);std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot open acoustic field cache");char magic[8];in.read(magic,8);if(!in||std::memcmp(magic,kMagic,8)!=0)throw std::runtime_error("invalid acoustic field cache magic/version");std::uint8_t h[16];in.read(reinterpret_cast<char*>(h),16);if(!in)throw std::runtime_error("truncated acoustic field cache header");std::vector<std::uint8_t> hb(h,h+16);Reader hr{hb};const auto size=hr.u64(),hash=hr.u64();if(size>kMaxPayload)throw std::runtime_error("acoustic field cache payload exceeds limit");std::vector<std::uint8_t> payload(static_cast<std::size_t>(size));in.read(reinterpret_cast<char*>(payload.data()),static_cast<std::streamsize>(payload.size()));if(!in||in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("truncated cache or unexpected trailing data");if(checksum(payload)!=hash)throw std::runtime_error("acoustic field cache checksum mismatch");
 Reader r{payload};if(r.u32()!=kFormatVersion)throw std::runtime_error("unsupported acoustic field cache format version");AcousticFieldCache cache;cache.identity={r.str(),r.str(),r.str()};if(cache.identity.geometry_fingerprint!=expected.geometry_fingerprint||cache.identity.material_fingerprint!=expected.material_fingerprint||cache.identity.solver_fingerprint!=expected.solver_fingerprint)throw std::runtime_error("acoustic field cache identity mismatch (geometry/material/solver)");const auto modes=r.u32();if(modes==0||modes>kMaxModes)throw std::runtime_error("invalid mode count in acoustic cache");cache.modes.reserve(modes);
 for(std::uint32_t m=0;m<modes;++m){CachedModalRadiationField entry;entry.frequency_hz=r.real();const auto nodes=r.u32();if(entry.frequency_hz<=0||nodes<3||nodes>kMaxNodes)throw std::runtime_error("invalid mode/polar dimensions in acoustic cache");auto& pf=entry.directional_field;pf.polar_angles_rad.reserve(nodes);pf.fields.reserve(nodes);
  for(std::uint32_t n=0;n<nodes;++n){const double polar=r.real(),radius=r.real(),obsPolar=r.real();const auto count=r.u32();if(polar<0||polar>kPi||radius<=0||obsPolar<0||obsPolar>kPi||count==0||count>kMaxCoefficients)throw std::runtime_error("invalid field dimensions in acoustic cache");AcousticDirectionalField f;f.observation_radius_m=radius;f.observation_polar_rad=obsPolar;f.orders.reserve(count);f.coefficients.reserve(count);
   for(std::uint32_t k=0;k<count;++k){const auto bits=r.u32();std::int32_t order=0;std::memcpy(&order,&bits,4);const double re=r.real(),im=r.real();f.orders.push_back(order);f.coefficients.emplace_back(re,im);}if(!f.valid())throw std::runtime_error("invalid harmonic field in acoustic cache");pf.polar_angles_rad.push_back(polar);pf.fields.push_back(std::move(f));}
  if(!pf.valid()) throw std::runtime_error("invalid polar field in acoustic cache");
  cache.modes.push_back(std::move(entry));}
 if(r.p!=payload.size()||!cache.valid()) throw std::runtime_error("invalid or trailing acoustic cache data");
 return cache;
}

PreparedStereoRadiation prepareStereoRadiationFromCache(const AcousticFieldCache& cache,double distance,double polar,double azimuth,double spacing,double maxYaw){
 if(!cache.valid()||!(std::isfinite(distance)&&distance>0.25)||!std::isfinite(polar)||polar<0||polar>kPi||!std::isfinite(azimuth)||!(std::isfinite(spacing)&&spacing>0&&spacing<0.5)||!(std::isfinite(maxYaw)&&maxYaw>=0&&maxYaw<=0.5*kPi))throw std::invalid_argument("invalid cached stereo receiver configuration");
 PreparedStereoRadiation prepared;prepared.centreDistanceM_=distance;prepared.centrePolarRad_=polar;prepared.centreAzimuthRad_=azimuth;prepared.earSpacingM_=spacing;prepared.maximumAbsYawRad_=maxYaw;prepared.fields_.reserve(cache.modes.size());prepared.mode_frequencies_hz_.reserve(cache.modes.size());
 for(const auto& mode:cache.modes){prepared.fields_.push_back(mode.directional_field);prepared.mode_frequencies_hz_.push_back(mode.frequency_hz);}return prepared;
}

} // namespace BellOptimizer::Acoustics
