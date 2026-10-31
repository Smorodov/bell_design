#include "cylindrical_q4_element.hpp"

#include <cmath>
#include <stdexcept>

namespace BellOptimizer::FEM {
namespace {
constexpr double PI = 3.141592653589793238462643383279502884;
using namespace BellOptimizer::Math;

struct BasisPoint { double value, dphi; };
BasisPoint angular_basis(int harmonic, bool sine, double phi) {
    if (harmonic == 0) return {1.0, 0.0};
    const double h = static_cast<double>(harmonic);
    return sine ? BasisPoint{std::sin(h*phi), h*std::cos(h*phi)}
                : BasisPoint{std::cos(h*phi), -h*std::sin(h*phi)};
}
}

CylindricalQ4ElementMatrices assemble_cylindrical_q4(
    const std::array<MeridionalNode, 4>& nodes, double E, double nu, double rho,
    const std::vector<int>& harmonics, unsigned angular_samples, double ellipticity) {
    if (!(ellipticity > 0.0 && ellipticity <= 1.0) || !std::isfinite(ellipticity))
        throw std::invalid_argument("cylindrical Q4: ellipticity must be finite in (0,1]");
    if (!(rho > 0.0) || !std::isfinite(rho))
        throw std::invalid_argument("cylindrical Q4: density must be finite and positive");
    if (angular_samples < 16 || angular_samples % 2 != 0)
        throw std::invalid_argument("cylindrical Q4: angular sample count must be even and >=16");
    if (harmonics.empty() || harmonics.front() != 0)
        throw std::invalid_argument("cylindrical Q4: harmonics must start with zero");
    for (std::size_t i=0; i<harmonics.size(); ++i) {
        if (harmonics[i] < 0 || (i && harmonics[i] <= harmonics[i-1]))
            throw std::invalid_argument("cylindrical Q4: harmonics must be increasing, unique and nonnegative");
    }
    const auto C = isotropic_elasticity_3d(E, nu);
    CylindricalQ4ElementMatrices out;
    for (int h : harmonics) {
        out.harmonics.push_back(h); out.sine_basis.push_back(false);
        if (h) { out.harmonics.push_back(h); out.sine_basis.push_back(true); }
    }
    const std::size_t nb = out.harmonics.size();
    const std::size_t ndof = 12u * nb;
    out.stiffness.assign(ndof, std::vector<double>(ndof, 0.0));
    out.mass.assign(ndof, std::vector<double>(ndof, 0.0));

    const double g = 1.0/std::sqrt(3.0);
    const double gp[2] = {-g, g};
    const double dphi = 2.0*PI/static_cast<double>(angular_samples);
    const double inv_q2 = 1.0 / (ellipticity * ellipticity);
    for (double xi : gp) for (double eta : gp) {
        const double N[4] = {0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta),
                             0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta)};
        const double dNdxi[4] = {-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)};
        const double dNdeta[4] = {-0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi)};
        double J00=0,J01=0,J10=0,J11=0,r=0,z=0;
        for (int a=0;a<4;++a) {
            J00 += dNdxi[a]*nodes[a].r; J01 += dNdeta[a]*nodes[a].r;
            J10 += dNdxi[a]*nodes[a].z; J11 += dNdeta[a]*nodes[a].z;
            r += N[a]*nodes[a].r; z += N[a]*nodes[a].z;
        }
        (void)z;
        const double det = J00*J11-J01*J10;
        if (!(r > 0.0) || !(det > 0.0) || !std::isfinite(det))
            throw std::invalid_argument("cylindrical Q4: invalid orientation or quadrature point at r<=0");
        double dNdr[4], dNdz[4];
        for (int a=0;a<4;++a) {
            dNdr[a] = ( J11*dNdxi[a]-J10*dNdeta[a])/det;
            dNdz[a] = (-J01*dNdxi[a]+J00*dNdeta[a])/det;
        }
        for (unsigned ip=0; ip<angular_samples; ++ip) {
            const double phi=dphi*static_cast<double>(ip);
            double alpha = 1.0, alpha_log_derivative = 0.0;
            if (ellipticity != 1.0) {
                const double cphi = std::cos(phi), sphi = std::sin(phi);
                const double metric = cphi*cphi + sphi*sphi*inv_q2;
                alpha = 1.0 / std::sqrt(metric);
                alpha_log_derivative = -sphi*cphi*(inv_q2 - 1.0) / metric;
            }
            const double physical_r = alpha * r;
            const double physical_det = alpha * det;
            std::vector<BasisPoint> av; av.reserve(nb);
            for (int h : harmonics) {
                av.push_back(angular_basis(h,false,phi));
                if (h) av.push_back(angular_basis(h,true,phi));
            }
            std::vector<CylindricalStrain> B(ndof);
            std::vector<double> shape(ndof,0.0);
            for (int a=0;a<4;++a) for (int c=0;c<3;++c) for (std::size_t b=0;b<nb;++b) {
                const std::size_t dof=(static_cast<std::size_t>(a)*3u+static_cast<std::size_t>(c))*nb+b;
                const double v=N[a]*av[b].value;
                shape[dof]=v;
                CylindricalDisplacementJet jet; jet.r=r;
                const double radial=N[a]*av[b].value;
                // The physical azimuthal derivative is taken at fixed physical
                // radius. Since r=alpha(phi)*r_ref, dr_ref/dphi|r =
                // -(alpha'/alpha)*r_ref; this geometry term is essential for
                // a consistent non-axisymmetric strain-displacement operator.
                const double azimuthal=N[a]*av[b].dphi -
                    alpha_log_derivative*r*dNdr[a]*av[b].value;
                const double dr=(dNdr[a]/alpha)*av[b].value, dz=dNdz[a]*av[b].value;
                if (c==0) { jet.ur=radial; jet.dur_dr=dr; jet.dur_dphi=azimuthal; jet.dur_dz=dz; }
                else if (c==1) { jet.uphi=radial; jet.duphi_dr=dr; jet.duphi_dphi=azimuthal; jet.duphi_dz=dz; }
                else { jet.uz=radial; jet.duz_dr=dr; jet.duz_dphi=azimuthal; jet.duz_dz=dz; }
                jet.r=physical_r;
                B[dof]=cylindrical_small_strain(jet);
            }
            const double weight=physical_r*physical_det*dphi;
            for (std::size_t i=0;i<ndof;++i) for (std::size_t j=0;j<ndof;++j) {
                double kb=0.0;
                for (std::size_t p=0;p<6;++p) for (std::size_t q=0;q<6;++q)
                    kb += B[i][p]*C[p][q]*B[j][q];
                out.stiffness[i][j] += kb*weight;
                const std::size_t ci=(i/nb)%3u, cj=(j/nb)%3u;
                if (ci==cj) out.mass[i][j] += rho*shape[i]*shape[j]*weight;
            }
        }
    }
    return out;
}
} // namespace BellOptimizer::FEM
