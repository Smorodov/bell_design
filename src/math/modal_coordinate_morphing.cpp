#include "modal_coordinate_morphing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BellOptimizer::Math {
namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double SQRT3_INV = 0.57735026918962576451;

struct FieldSample {
    double ur = 0.0, uz = 0.0;
    double urr = 0.0, urz = 0.0;
    double uzr = 0.0, uzz = 0.0;
};

FieldSample interpolate_q4(const std::array<Node, 4>& nd,
                           const std::array<double, 4>& ur_nodes,
                           const std::array<double, 4>& uz_nodes,
                           double xi, double eta) {
    const double N[4] = {
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
    const double dNdxi[4] = {
        -0.25 * (1.0 - eta), 0.25 * (1.0 - eta),
        0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
    const double dNdeta[4] = {
        -0.25 * (1.0 - xi), -0.25 * (1.0 + xi),
        0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};

    double J00 = 0.0, J01 = 0.0, J10 = 0.0, J11 = 0.0;
    for (int a = 0; a < 4; ++a) {
        J00 += dNdxi[a] * nd[a].r;
        J01 += dNdxi[a] * nd[a].z;
        J10 += dNdeta[a] * nd[a].r;
        J11 += dNdeta[a] * nd[a].z;
    }
    const double det = J00 * J11 - J01 * J10;
    if (!(det > 0.0) || !std::isfinite(det))
        throw std::runtime_error("coordinate morphing: invalid Q4 Jacobian");

    const double inv00 = J11 / det;
    const double inv01 = -J01 / det;
    const double inv10 = -J10 / det;
    const double inv11 = J00 / det;
    double dNdr[4], dNdz[4];
    for (int a = 0; a < 4; ++a) {
        dNdr[a] = dNdxi[a] * inv00 + dNdeta[a] * inv01;
        dNdz[a] = dNdxi[a] * inv10 + dNdeta[a] * inv11;
    }

    FieldSample s;
    for (int a = 0; a < 4; ++a) {
        s.ur += N[a] * ur_nodes[a];
        s.uz += N[a] * uz_nodes[a];
        s.urr += dNdr[a] * ur_nodes[a];
        s.urz += dNdz[a] * ur_nodes[a];
        s.uzr += dNdr[a] * uz_nodes[a];
        s.uzz += dNdz[a] * uz_nodes[a];
    }

    double r = 0.0;
    for (int a = 0; a < 4; ++a) r += N[a] * nd[a].r;
    if (!(r > 0.0) || !std::isfinite(r))
        throw std::runtime_error("coordinate morphing: non-positive integration radius");
    return s;
}

struct EnergyPair {
    double mass = 0.0;
    double strain = 0.0;
};

EnergyPair integrate_mode(const BellMeshData& mesh,
                           const BellMaterial& material,
                           const ModalMode& mode,
                           const CoordinateScale3D& scale) {
    if (mode.nodal_radial.size() != mesh.nodes.size() ||
        mode.nodal_axial.size() != mesh.nodes.size())
        throw std::invalid_argument(
            "coordinate morphing: mode nodal fields are required; solve with capture_mode_shapes=true");
    if (!(material.E > 0.0) || !(material.rho > 0.0) ||
        !std::isfinite(material.E) || !std::isfinite(material.rho) ||
        !std::isfinite(material.nu) || material.nu <= -1.0 || material.nu >= 0.5)
        throw std::invalid_argument("coordinate morphing: invalid material");

    // 3-D isotropic elasticity matrix in Voigt order xx, yy, zz, xy, yz, xz,
    // with engineering shear strains. The base FEM uses the same Lamé law.
    const double lambda = material.E * material.nu /
        ((1.0 + material.nu) * (1.0 - 2.0 * material.nu));
    const double mu = material.E / (2.0 * (1.0 + material.nu));
    double D[6][6]{};
    D[0][0] = D[1][1] = D[2][2] = lambda + 2.0 * mu;
    D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = lambda;
    D[3][3] = D[4][4] = D[5][5] = mu;

    EnergyPair out;
    // 2x2 Q4 quadrature in the meridional element and 8 equally spaced azimuth
    // points. The transformed quadratic strain energy is exactly integrated in
    // azimuth for this affine transform; eight points provide a transparent
    // numerical implementation and a useful convergence check.
    constexpr int phi_count = 8;
    std::vector<int> id_to_index(mesh.nodes.size() + 1u, -1);
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const int id = mesh.nodes[i].id;
        if (id <= 0 || static_cast<size_t>(id) >= id_to_index.size() || id_to_index[static_cast<size_t>(id)] >= 0)
            throw std::invalid_argument("coordinate morphing: invalid node IDs");
        id_to_index[static_cast<size_t>(id)] = static_cast<int>(i);
    }

    for (const auto& e : mesh.elements) {
        const int element_ids[4] = {e.n1, e.n2, e.n3, e.n4};
        int ids[4]{};
        for (int a = 0; a < 4; ++a) {
            if (element_ids[a] <= 0 || static_cast<size_t>(element_ids[a]) >= id_to_index.size())
                throw std::invalid_argument("coordinate morphing: invalid element node");
            ids[a] = id_to_index[static_cast<size_t>(element_ids[a])];
            if (ids[a] < 0) throw std::invalid_argument("coordinate morphing: unknown element node");
        }

        std::array<Node, 4> nd;
        std::array<double, 4> urn{}, uzn{};
        for (int a = 0; a < 4; ++a) {
            nd[a] = mesh.nodes[static_cast<size_t>(ids[a])];
            urn[a] = mode.nodal_radial[static_cast<size_t>(ids[a])];
            uzn[a] = mode.nodal_axial[static_cast<size_t>(ids[a])];
        }

        for (double xi : {-SQRT3_INV, SQRT3_INV}) {
            for (double eta : {-SQRT3_INV, SQRT3_INV}) {
                // Rebuild the local geometry Jacobian for the integration weight.
                const double dNdxi[4] = {
                    -0.25 * (1.0 - eta), 0.25 * (1.0 - eta),
                    0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
                const double dNdeta[4] = {
                    -0.25 * (1.0 - xi), -0.25 * (1.0 + xi),
                    0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};
                const double N[4] = {
                    0.25 * (1.0 - xi) * (1.0 - eta),
                    0.25 * (1.0 + xi) * (1.0 - eta),
                    0.25 * (1.0 + xi) * (1.0 + eta),
                    0.25 * (1.0 - xi) * (1.0 + eta)};
                double J00 = 0.0, J01 = 0.0, J10 = 0.0, J11 = 0.0;
                double r = 0.0;
                for (int a = 0; a < 4; ++a) {
                    J00 += dNdxi[a] * nd[a].r;
                    J01 += dNdxi[a] * nd[a].z;
                    J10 += dNdeta[a] * nd[a].r;
                    J11 += dNdeta[a] * nd[a].z;
                    r += N[a] * nd[a].r;
                }
                const double det = J00 * J11 - J01 * J10;
                if (!(det > 0.0) || !(r > 0.0) || !std::isfinite(det) || !std::isfinite(r))
                    throw std::runtime_error("coordinate morphing: invalid integration point");

                const FieldSample f = interpolate_q4(nd, urn, uzn, xi, eta);
                const double ur_over_r = f.ur / r;

                for (int ip = 0; ip < phi_count; ++ip) {
                    const double phi = 2.0 * PI * static_cast<double>(ip) / phi_count;
                    const double c = std::cos(phi), s = std::sin(phi);

                    // Base cylindrical displacement gradient in Cartesian axes.
                    double G[3][3]{};
                    G[0][0] = f.urr * c * c + ur_over_r * s * s;
                    G[0][1] = (f.urr - ur_over_r) * s * c;
                    G[1][0] = G[0][1];
                    G[1][1] = f.urr * s * s + ur_over_r * c * c;
                    G[0][2] = f.urz * c;
                    G[1][2] = f.urz * s;
                    G[2][0] = f.uzr * c;
                    G[2][1] = f.uzr * s;
                    G[2][2] = f.uzz;

                    // Affine coordinate/displacement transformation: G' = S G S^-1.
                    const double sv[3] = {scale.x, scale.y, scale.z};
                    double Gp[3][3]{};
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            Gp[i][j] = sv[i] * G[i][j] / sv[j];

                    const double exx = Gp[0][0];
                    const double eyy = Gp[1][1];
                    const double ezz = Gp[2][2];
                    const double exy = Gp[0][1] + Gp[1][0];
                    const double eyz = Gp[1][2] + Gp[2][1];
                    const double exz = Gp[0][2] + Gp[2][0];
                    const double strain[6] = {exx, eyy, ezz, exy, eyz, exz};

                    double stress[6]{};
                    for (int i = 0; i < 6; ++i)
                        for (int j = 0; j < 6; ++j)
                            stress[i] += D[i][j] * strain[j];
                    double density_energy = 0.0;
                    for (int i = 0; i < 6; ++i) density_energy += 0.5 * strain[i] * stress[i];

                    const double ux = scale.x * f.ur * c;
                    const double uy = scale.y * f.ur * s;
                    const double uz = scale.z * f.uz;
                    const double displacement_sq = ux * ux + uy * uy + uz * uz;

                    const double weight = (2.0 * PI / phi_count) * r * det * scale.volume_scale();
                    out.strain += density_energy * weight;
                    out.mass += material.rho * displacement_sq * weight;
                }
            }
        }
    }
    return out;
}

} // namespace

MechanicalMorphingResult evaluate_mechanical_morphing(
    const BellMeshData& mesh,
    const BellMaterial& material,
    const ModalMode& mode,
    const CoordinateScale3D& scale) {
    scale.validate();
    const EnergyPair base = integrate_mode(mesh, material, mode, CoordinateScale3D{});
    const EnergyPair morphed = integrate_mode(mesh, material, mode, scale);
    if (!(base.mass > 0.0) || !(base.strain > 0.0) ||
        !(morphed.mass > 0.0) || !(morphed.strain > 0.0))
        throw std::runtime_error("coordinate morphing: non-positive modal energy");

    MechanicalMorphingResult out;
    out.base_mass = base.mass;
    out.base_strain_energy = base.strain;
    out.morphed_mass = morphed.mass;
    out.morphed_strain_energy = morphed.strain;
    out.base_frequency_hz = std::sqrt(2.0 * base.strain / base.mass) / (2.0 * PI);
    out.morphed_frequency_hz = std::sqrt(2.0 * morphed.strain / morphed.mass) / (2.0 * PI);
    out.relative_mass_change = morphed.mass / base.mass - 1.0;
    out.relative_energy_change = morphed.strain / base.strain - 1.0;
    out.relative_frequency_change = out.morphed_frequency_hz / out.base_frequency_hz - 1.0;
    return out;
}

} // namespace BellOptimizer::Math
