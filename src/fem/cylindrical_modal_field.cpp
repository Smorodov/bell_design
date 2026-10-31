#include "cylindrical_modal_field.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace BellOptimizer::FEM {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss_x[3] = {-0.77459666924148337704, 0.0, 0.77459666924148337704};
constexpr double gauss_w[3] = {0.55555555555555555556, 0.88888888888888888889, 0.55555555555555555556};

struct ShapePoint {
    std::array<double, 4> n{};
    double dr_dxi = 0.0, dr_deta = 0.0;
    double dz_dxi = 0.0, dz_deta = 0.0;
    double det_j = 0.0;
};

ShapePoint shape_at(const BellMeshData& mesh, const Element& e, double xi, double eta) {
    const int ids[4] = {e.n1, e.n2, e.n3, e.n4};
    ShapePoint out;
    const double dx[4] = {-0.25*(1.0-eta), 0.25*(1.0-eta), 0.25*(1.0+eta), -0.25*(1.0+eta)};
    const double de[4] = {-0.25*(1.0-xi), -0.25*(1.0+xi), 0.25*(1.0+xi), 0.25*(1.0-xi)};
    double dr_dxi=0.0, dr_deta=0.0, dz_dxi=0.0, dz_deta=0.0;
    for (int i=0;i<4;++i) {
        if (ids[i] <= 0 || static_cast<std::size_t>(ids[i]) > mesh.nodes.size())
            throw std::out_of_range("modal field: element node index out of range");
        const auto& p = mesh.nodes[static_cast<std::size_t>(ids[i]-1)];
        out.n[static_cast<std::size_t>(i)] = 0.25 *
            (1.0 + (i==1 || i==2 ? xi : -xi)) *
            (1.0 + (i>=2 ? eta : -eta));
        dr_dxi += dx[i]*p.r; dr_deta += de[i]*p.r;
        dz_dxi += dx[i]*p.z; dz_deta += de[i]*p.z;
    }
    out.dr_dxi=dr_dxi; out.dr_deta=dr_deta;
    out.dz_dxi=dz_dxi; out.dz_deta=dz_deta;
    out.det_j = dr_dxi*dz_deta - dr_deta*dz_dxi;
    return out;
}

std::pair<double,double> map_point(const BellMeshData& mesh, const Element& e,
                                   double xi, double eta) {
    const auto s = shape_at(mesh,e,xi,eta);
    const int ids[4] = {e.n1,e.n2,e.n3,e.n4};
    double r=0.0,z=0.0;
    for(int i=0;i<4;++i) {
        const auto& p=mesh.nodes[static_cast<std::size_t>(ids[i]-1)];
        r += s.n[static_cast<std::size_t>(i)]*p.r;
        z += s.n[static_cast<std::size_t>(i)]*p.z;
    }
    return {r,z};
}

struct ElementBox { double r0,r1,z0,z1; };

class ElementLocator {
public:
    explicit ElementLocator(const BellMeshData& mesh) : mesh_(mesh) {
        if(mesh.nodes.empty() || mesh.elements.empty()) throw std::invalid_argument("modal field: empty mesh");
        r0_=z0_=std::numeric_limits<double>::infinity();
        r1_=z1_=-std::numeric_limits<double>::infinity();
        for(const auto& p:mesh.nodes) { r0_=std::min(r0_,p.r); r1_=std::max(r1_,p.r); z0_=std::min(z0_,p.z); z1_=std::max(z1_,p.z); }
        const std::size_t side=std::max<std::size_t>(8,std::min<std::size_t>(128,
            static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(mesh.elements.size()))*2.0))));
        nr_=side; nz_=side;
        if (nr_ != 0 && nz_ > std::numeric_limits<std::size_t>::max() / nr_)
            throw std::length_error("modal field: locator grid size overflow");
        bins_.resize(nr_*nz_);
        if (mesh.elements.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("modal field: element count exceeds locator index range");
        boxes_.reserve(mesh.elements.size());
        for(std::size_t k=0;k<mesh.elements.size();++k) {
            const auto& e=mesh.elements[k]; const int ids[4]={e.n1,e.n2,e.n3,e.n4};
            ElementBox b{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
            for(int id:ids) {
                if(id<=0 || static_cast<std::size_t>(id)>mesh.nodes.size()) throw std::out_of_range("modal field: invalid mesh connectivity");
                const auto& p=mesh.nodes[static_cast<std::size_t>(id-1)];
                b.r0=std::min(b.r0,p.r); b.r1=std::max(b.r1,p.r); b.z0=std::min(b.z0,p.z); b.z1=std::max(b.z1,p.z);
            }
            boxes_.push_back(b);
            const auto ir0=bin_r(b.r0), ir1=bin_r(b.r1), iz0=bin_z(b.z0), iz1=bin_z(b.z1);
            for(std::size_t iz=iz0;iz<=iz1;++iz) for(std::size_t ir=ir0;ir<=ir1;++ir)
                bins_[iz*nr_+ir].push_back(static_cast<int>(k));
        }
    }

    bool locate(double r,double z,int& element,double& xi,double& eta,ShapePoint& shape) const {
        const double eps=1e-10*std::max({1.0,std::abs(r1_-r0_),std::abs(z1_-z0_)});
        if(r<r0_-eps || r>r1_+eps || z<z0_-eps || z>z1_+eps) return false;
        const auto& candidates=bins_[bin_z(z)*nr_+bin_r(r)];
        for(int index:candidates) {
            const auto& b=boxes_[static_cast<std::size_t>(index)];
            if(r<b.r0-eps || r>b.r1+eps || z<b.z0-eps || z>b.z1+eps) continue;
            const auto& e=mesh_.elements[static_cast<std::size_t>(index)];
            double u=0.0,v=0.0;
            for(int it=0;it<16;++it) {
                const auto p=map_point(mesh_,e,u,v);
                const double fr=p.first-r, fz=p.second-z;
                if(std::max(std::abs(fr),std::abs(fz))<=eps) break;
                const auto jacobian=shape_at(mesh_,e,u,v);
                const double det=jacobian.det_j;
                if(!(det>1e-18) || !std::isfinite(det)) break;
                const double du=(jacobian.dz_deta*fr-jacobian.dr_deta*fz)/det;
                const double dv=(-jacobian.dz_dxi*fr+jacobian.dr_dxi*fz)/det;
                u-=du; v-=dv;
                if(!std::isfinite(u)||!std::isfinite(v)) break;
            }
            if(u>=-1.0-1e-8 && u<=1.0+1e-8 && v>=-1.0-1e-8 && v<=1.0+1e-8) {
                const auto p=map_point(mesh_,e,u,v);
                if(std::hypot(p.first-r,p.second-z)<=eps*10.0) {
                    element=index; xi=std::max(-1.0,std::min(1.0,u)); eta=std::max(-1.0,std::min(1.0,v));
                    shape=shape_at(mesh_,e,xi,eta);
                    if(!(shape.det_j>0.0) || !std::isfinite(shape.det_j)) return false;
                    return true;
                }
            }
        }
        return false;
    }
private:
    std::size_t bin_r(double r) const { return coordinate_bin(r,r0_,r1_,nr_); }
    std::size_t bin_z(double z) const { return coordinate_bin(z,z0_,z1_,nz_); }
    static std::size_t coordinate_bin(double x,double lo,double hi,std::size_t n) {
        if(!(hi>lo)) return 0;
        const double t=std::max(0.0,std::min(1.0,(x-lo)/(hi-lo)));
        return std::min(n-1,static_cast<std::size_t>(t*static_cast<double>(n)));
    }
    const BellMeshData& mesh_; double r0_=0,r1_=0,z0_=0,z1_=0;
    std::size_t nr_=0,nz_=0; std::vector<ElementBox> boxes_; std::vector<std::vector<int>> bins_;
};

std::size_t basis_count(const std::vector<int>& h) {
    if(h.empty() || h.front()!=0) throw std::invalid_argument("modal field: harmonics must start at zero");
    std::size_t n=0; int previous=-1;
    for(int v:h) { if(v<0 || v<=previous) throw std::invalid_argument("modal field: harmonics must be strictly increasing and nonnegative"); n+=(v==0?1u:2u); previous=v; }
    return n;
}

void validate_field(const CylindricalModalField& f) {
    if(!f.mesh || !f.harmonics || !f.mode_shape) throw std::invalid_argument("modal field: null view member");
    const auto nb=basis_count(*f.harmonics);
    if (f.mesh->nodes.size() > std::numeric_limits<std::size_t>::max() / 3u ||
        nb > std::numeric_limits<std::size_t>::max() / (f.mesh->nodes.size() * 3u))
        throw std::length_error("modal field: modal vector size overflow");
    const std::size_t expected = f.mesh->nodes.size() * 3u * nb;
    if(f.mesh->nodes.empty() || f.mesh->elements.empty() || f.mode_shape->size()!=expected)
        throw std::invalid_argument("modal field: inconsistent mesh, harmonics or modal vector size");
    if(!(f.ellipticity>0.0 && f.ellipticity<=1.0) || !std::isfinite(f.ellipticity) ||
       !(f.density_kg_m3>0.0) || !std::isfinite(f.density_kg_m3))
        throw std::invalid_argument("modal field: invalid ellipticity or density");
}

std::array<double,3> evaluate(const CylindricalModalField& f,const Element& e,const ShapePoint& s,
                              const std::vector<double>& angular) {
    const auto& v=*f.mode_shape;
    const std::size_t nb=angular.size(); std::array<double,3> u{{0.0,0.0,0.0}};
    const int ids[4]={e.n1,e.n2,e.n3,e.n4};
    for(int a=0;a<4;++a) {
        if (ids[a] <= 0 || static_cast<std::size_t>(ids[a]) > f.mesh->nodes.size())
            throw std::out_of_range("modal field: element node index out of range during evaluation");
        const std::size_t node=static_cast<std::size_t>(ids[a]-1);
        const double n=s.n[static_cast<std::size_t>(a)];
        for(std::size_t c=0;c<3;++c) {
            double value=0.0;
            for(std::size_t b=0;b<nb;++b) value+=v[(node*3u+c)*nb+b]*angular[b];
            u[c]+=n*value;
        }
    }
    return u;
}

struct Integral { double mass=0.0; double volume=0.0; };

Integral integrate_direction(const CylindricalModalField& integration,
                              const CylindricalModalField& evaluation,
                              const ElementLocator* locator,unsigned angular_samples) {
    const auto& mesh=*integration.mesh; const auto& eval_mesh=*evaluation.mesh;
    const double dphi=2.0*pi/static_cast<double>(angular_samples);
    const double density=0.5*(integration.density_kg_m3+evaluation.density_kg_m3);
    const double invq2=1.0/(integration.ellipticity*integration.ellipticity);
    std::vector<std::vector<double>> angular_values(angular_samples), evaluation_angular_values(angular_samples);
    for(unsigned ip=0;ip<angular_samples;++ip) {
        const double phi=dphi*static_cast<double>(ip);
        for(int harmonic:*integration.harmonics) {
            angular_values[ip].push_back(harmonic==0?1.0:std::cos(static_cast<double>(harmonic)*phi));
            if(harmonic!=0) angular_values[ip].push_back(std::sin(static_cast<double>(harmonic)*phi));
        }
        for(int harmonic:*evaluation.harmonics) {
            evaluation_angular_values[ip].push_back(harmonic==0?1.0:std::cos(static_cast<double>(harmonic)*phi));
            if(harmonic!=0) evaluation_angular_values[ip].push_back(std::sin(static_cast<double>(harmonic)*phi));
        }
    }
    Integral result;
    for(const auto& e:mesh.elements) {
        for(int ir=0;ir<3;++ir) for(int iz=0;iz<3;++iz) {
            const double xi=gauss_x[ir], eta=gauss_x[iz];
            const auto shape=shape_at(mesh,e,xi,eta);
            if(!(shape.det_j>0.0) || !std::isfinite(shape.det_j)) throw std::invalid_argument("modal field: invalid Q4 orientation");
            const int ids[4]={e.n1,e.n2,e.n3,e.n4}; double r=0.0,z=0.0;
            for(int a=0;a<4;++a) { const auto& p=mesh.nodes[static_cast<std::size_t>(ids[a]-1)]; r+=shape.n[static_cast<std::size_t>(a)]*p.r; z+=shape.n[static_cast<std::size_t>(a)]*p.z; }
            if(!(r>0.0)) throw std::invalid_argument("modal field: nonpositive radius");
            int target_element=-1; double tx=0.0,ty=0.0; ShapePoint target_shape;
            if(locator && !locator->locate(r,z,target_element,tx,ty,target_shape)) continue;
            const Element& eval_element=locator?eval_mesh.elements[static_cast<std::size_t>(target_element)]:e;
            if(!locator) target_shape=shape;
            const double base_weight=gauss_w[ir]*gauss_w[iz]*shape.det_j*r;
            for(unsigned ip=0;ip<angular_samples;++ip) {
                const double phi=dphi*static_cast<double>(ip);
                double alpha=1.0;
                if(integration.ellipticity!=1.0) {
                    const double c=std::cos(phi),s=std::sin(phi);
                    alpha=1.0/std::sqrt(c*c+s*s*invq2);
                }
                const double weight=base_weight*alpha*alpha*dphi;
                const auto& angular=angular_values[ip];
                const auto ua=evaluate(integration,e,shape,angular);
                const auto ub=evaluate(evaluation,eval_element,target_shape,evaluation_angular_values[ip]);
                result.mass+=density*(ua[0]*ub[0]+ua[1]*ub[1]+ua[2]*ub[2])*weight;
                result.volume+=weight;
            }
        }
    }
    return result;
}

double inner_product(const CylindricalModalField& a,const CylindricalModalField& b,
                     unsigned angular_samples,double& coverage) {
    validate_field(a); validate_field(b);
    if(std::abs(a.ellipticity-b.ellipticity)>1e-10 ||
       std::abs(a.density_kg_m3-b.density_kg_m3)>1e-8*std::max(a.density_kg_m3,b.density_kg_m3))
        throw std::invalid_argument("modal field: fields must share ellipticity and density for a common mass metric");
    const int max_a=*std::max_element(a.harmonics->begin(),a.harmonics->end());
    const int max_b=*std::max_element(b.harmonics->begin(),b.harmonics->end());
    if(static_cast<long long>(max_a)+static_cast<long long>(max_b)>=static_cast<long long>(angular_samples))
        throw std::invalid_argument("modal field: angular quadrature aliases the product of Fourier bases");
    if(a.mesh==b.mesh) {
        const auto ab=integrate_direction(a,b,nullptr,angular_samples);
        coverage=1.0;
        return ab.mass;
    }
    ElementLocator locator_b(*b.mesh), locator_a(*a.mesh);
    const auto ab=integrate_direction(a,b,&locator_b,angular_samples);
    const auto ba=integrate_direction(b,a,&locator_a,angular_samples);
    const double va=integrate_direction(a,a,nullptr,angular_samples).volume;
    const double vb=integrate_direction(b,b,nullptr,angular_samples).volume;
    coverage=std::min(ab.volume/std::max(va,1e-300),ba.volume/std::max(vb,1e-300));
    return 0.5*(ab.mass+ba.mass);
}

double self_product(const CylindricalModalField& a,unsigned angular_samples) {
    return integrate_direction(a,a,nullptr,angular_samples).mass;
}

double smallest_singular_value_2x2(double a,double b,double c,double d) {
    const double tr=a*a+b*b+c*c+d*d, det=a*d-b*c;
    const double disc=std::max(0.0,tr*tr-4.0*det*det);
    const double lambda_max=0.5*(tr+std::sqrt(disc));
    const double lambda_min=lambda_max>0.0?(det*det/lambda_max):0.0;
    return std::sqrt(std::max(0.0,lambda_min));
}

} // namespace

CylindricalModalFieldOverlap compare_cylindrical_modal_fields(
    const CylindricalModalField& a,const CylindricalModalField& b,
    unsigned angular_samples,double minimum_common_volume_fraction) {
    validate_field(a); validate_field(b);
    if(angular_samples<16 || angular_samples%2!=0) throw std::invalid_argument("modal field: angular samples must be even and >=16");
    if(!(minimum_common_volume_fraction>0.0 && minimum_common_volume_fraction<=1.0)) throw std::invalid_argument("modal field: invalid common-volume threshold");
    if(std::abs(a.ellipticity-b.ellipticity)>1e-10 ||
       std::abs(a.density_kg_m3-b.density_kg_m3)>1e-8*std::max(a.density_kg_m3,b.density_kg_m3))
        throw std::invalid_argument("modal field: fields must share ellipticity and density for a common mass metric");
    double coverage=0.0; const double cross=inner_product(a,b,angular_samples,coverage);
    const double na=self_product(a,angular_samples), nb=self_product(b,angular_samples);
    CylindricalModalFieldOverlap out; out.common_volume_fraction=coverage; out.cross_mass=cross;
    out.mass_norm_a=na; out.mass_norm_b=nb;
    if(!(na>0.0) || !(nb>0.0) || !std::isfinite(cross) || coverage<minimum_common_volume_fraction) return out;
    out.normalized_overlap=std::max(0.0,std::min(1.0,std::abs(cross)/std::sqrt(na*nb)));
    out.valid=true; return out;
}

double cylindrical_modal_subspace_similarity(
    const std::vector<CylindricalModalField>& a,const std::vector<CylindricalModalField>& b,
    unsigned angular_samples,double minimum_common_volume_fraction,double* common_volume_fraction) {
    if(a.empty() || a.size()!=b.size() || a.size()>2) throw std::invalid_argument("modal field: subspace dimension must be 1 or 2 and equal");
    if(angular_samples<16 || angular_samples%2!=0) throw std::invalid_argument("modal field: angular samples must be even and >=16");
    if(!(minimum_common_volume_fraction>0.0 && minimum_common_volume_fraction<=1.0)) throw std::invalid_argument("modal field: invalid common-volume threshold");
    for(const auto& field:a) validate_field(field);
    for(const auto& field:b) validate_field(field);
    double min_coverage=1.0;
    const std::size_t n=a.size(); double ga[2][2]={{0,0},{0,0}},gb[2][2]={{0,0},{0,0}},c[2][2]={{0,0},{0,0}};
    for(std::size_t i=0;i<n;++i) {
        ga[i][i]=self_product(a[i],angular_samples); gb[i][i]=self_product(b[i],angular_samples);
        for(std::size_t j=i+1;j<n;++j) {
            double coverage_a=0.0,coverage_b=0.0;
            ga[i][j]=ga[j][i]=inner_product(a[i],a[j],angular_samples,coverage_a);
            gb[i][j]=gb[j][i]=inner_product(b[i],b[j],angular_samples,coverage_b);
            min_coverage=std::min(min_coverage,std::min(coverage_a,coverage_b));
        }
        for(std::size_t j=0;j<n;++j) {
            double coverage=0.0; c[i][j]=inner_product(a[i],b[j],angular_samples,coverage); min_coverage=std::min(min_coverage,coverage);
        }
    }
    if(common_volume_fraction) *common_volume_fraction=min_coverage;
    if(min_coverage<minimum_common_volume_fraction) return 0.0;
    auto orthonormalize=[](double g[2][2],double q[2][2])->bool {
        if(!(g[0][0]>0.0)) return false;
        const double n0=std::sqrt(g[0][0]); q[0][0]=1.0/n0;
        if(g[0][1]!=0.0) q[1][0]=0.0;
        if(g[1][1]<=0.0) return false;
        if(g[0][1]==0.0) { q[0][1]=0.0; q[1][1]=1.0/std::sqrt(g[1][1]); return true; }
        const double alpha=g[0][1]/g[0][0];
        const double residual=g[1][1]-g[0][1]*g[0][1]/g[0][0];
        if(!(residual>1e-14*g[1][1])) return false;
        q[0][1]=-alpha/std::sqrt(residual); q[1][1]=1.0/std::sqrt(residual); return true;
    };
    if(n==1) {
        if(!(ga[0][0]>0.0 && gb[0][0]>0.0)) return 0.0;
        return std::min(1.0,std::abs(c[0][0])/std::sqrt(ga[0][0]*gb[0][0]));
    }
    double qa[2][2]={{0,0},{0,0}},qb[2][2]={{0,0},{0,0}};
    if(!orthonormalize(ga,qa) || !orthonormalize(gb,qb)) return 0.0;
    double cc[2][2]={{0,0},{0,0}};
    for(int i=0;i<2;++i) for(int j=0;j<2;++j)
        for(int p=0;p<2;++p) for(int q=0;q<2;++q) cc[i][j]+=qa[p][i]*c[p][q]*qb[q][j];
    return std::min(1.0,smallest_singular_value_2x2(cc[0][0],cc[0][1],cc[1][0],cc[1][1]));
}

} // namespace BellOptimizer::FEM
