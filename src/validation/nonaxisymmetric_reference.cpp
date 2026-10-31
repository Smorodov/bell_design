#include "nonaxisymmetric_reference.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace BellOptimizer::Validation {
namespace {
constexpr double PI = 3.141592653589793238462643383279502884;

void fail(std::string* error, const std::string& message) {
    if (error) *error = message;
}

} // namespace

NonAxisymmetricReferenceMesh build_reference_mesh(const BellMeshData& mesh,
                                                  double q,
                                                  int azimuthal_divisions) {
    if (!(q > 0.0) || !std::isfinite(q))
        throw std::invalid_argument("reference mesh: q must be finite and positive");
    if (azimuthal_divisions < 4)
        throw std::invalid_argument("reference mesh: at least four azimuthal divisions are required");
    if (mesh.nodes.empty() || mesh.elements.empty())
        throw std::invalid_argument("reference mesh: empty axisymmetric mesh");

    NonAxisymmetricReferenceMesh out;
    const int nr = static_cast<int>(mesh.nodes.size());
    out.nodes.reserve(static_cast<size_t>(nr) * static_cast<size_t>(azimuthal_divisions));
    for (int j = 0; j < azimuthal_divisions; ++j) {
        const double phi = 2.0 * PI * static_cast<double>(j) /
                           static_cast<double>(azimuthal_divisions);
        const double c = std::cos(phi), s = std::sin(phi);
        for (const auto& node : mesh.nodes) {
            out.nodes.push_back({
                node.id + j * nr,
                node.r * c,
                q * node.r * s,
                node.z});
        }
    }

    out.elements.reserve(mesh.elements.size() * static_cast<size_t>(azimuthal_divisions));
    for (int j = 0; j < azimuthal_divisions; ++j) {
        const int next = (j + 1) % azimuthal_divisions;
        const int a = j * nr;
        const int b = next * nr;
        for (const auto& e : mesh.elements) {
            ReferenceHex8 h{};
            // Bottom face follows the meridional Q4 order. The upper face is
            // repeated at the next azimuth with the same orientation.
            // C3D8 local orientation must have a positive Jacobian. The
            // meridional Q4 orientation is clockwise when viewed from the
            // positive azimuthal side, so reverse each face consistently.
            h.n[0] = e.n2 + a;
            h.n[1] = e.n1 + a;
            h.n[2] = e.n4 + a;
            h.n[3] = e.n3 + a;
            h.n[4] = e.n2 + b;
            h.n[5] = e.n1 + b;
            h.n[6] = e.n4 + b;
            h.n[7] = e.n3 + b;
            out.elements.push_back(h);
        }
    }

    // The affine y scaling changes volume by q exactly.
    // The axisymmetric volume is recovered from the meridional Q4 mesh by
    // polygon integration; this is only a topology/geometry audit quantity.
    double v2d = 0.0;
    for (const auto& e : mesh.elements) {
        const Node* n[4]{};
        for (int a = 0; a < 4; ++a) {
            const int id = (&e.n1)[a];
            if (id <= 0 || static_cast<size_t>(id) > mesh.nodes.size())
                throw std::invalid_argument("reference mesh: invalid element node id");
            n[a] = &mesh.nodes[static_cast<size_t>(id - 1)];
        }
        // 2x2 Gauss integration of 2*pi*r over the Q4 meridional area.
        constexpr double g = 0.57735026918962576451;
        for (double xi : {-g, g}) for (double eta : {-g, g}) {
            const double N[4] = {
                .25*(1-xi)*(1-eta), .25*(1+xi)*(1-eta),
                .25*(1+xi)*(1+eta), .25*(1-xi)*(1+eta)};
            const double dx[4] = {
                -.25*(1-eta), .25*(1-eta), .25*(1+eta), -.25*(1+eta)};
            const double de[4] = {
                -.25*(1-xi), -.25*(1+xi), .25*(1+xi), .25*(1-xi)};
            double j00=0,j01=0,j10=0,j11=0,r=0;
            for (int a=0;a<4;++a) {
                j00 += dx[a]*n[a]->r; j01 += dx[a]*n[a]->z;
                j10 += de[a]*n[a]->r; j11 += de[a]*n[a]->z;
                r += N[a]*n[a]->r;
            }
            const double det=j00*j11-j01*j10;
            if (!(det>0.0) || !(r>=0.0) || !std::isfinite(det))
                throw std::invalid_argument("reference mesh: invalid meridional element");
            v2d += 2.0*PI*r*det;
        }
    }
    out.computational_volume = v2d;
    out.physical_volume = q * v2d;
    return out;
}

bool write_calculix_c3d8_frequency_input(const NonAxisymmetricReferenceMesh& mesh,
                                         const BellMaterial& material,
                                         const std::string& filename,
                                         int requested_modes,
                                         std::string* error) {
    if (mesh.nodes.empty() || mesh.elements.empty()) {
        fail(error, "reference input: empty mesh"); return false;
    }
    if (!(material.E > 0.0) || !(material.rho > 0.0) ||
        !(material.nu > -1.0 && material.nu < 0.5)) {
        fail(error, "reference input: invalid material"); return false;
    }
    if (requested_modes < 1) {
        fail(error, "reference input: requested_modes must be positive"); return false;
    }
    std::ofstream f(filename);
    if (!f) { fail(error, "reference input: cannot open output file"); return false; }
    f << "*HEADING\n"
      << "BellOptimizer validation-only non-axisymmetric C3D8 reference\n"
      << "*MATERIAL, NAME=BRONZE\n*ELASTIC\n"
      << material.E << ", " << material.nu << "\n*DENSITY\n"
      << material.rho << "\n";
    f << "*NODE, NSET=NALL\n";
    for (const auto& n : mesh.nodes) f << n.id << ", " << n.x << ", " << n.y << ", " << n.z << "\n";
    f << "*ELEMENT, TYPE=C3D8, ELSET=EALL\n";
    for (size_t i=0;i<mesh.elements.size();++i) {
        const auto& e=mesh.elements[i];
        f << (i+1);
        for (int k=0;k<8;++k) f << ", " << e.n[k];
        f << '\n';
    }
    f << "*SOLID SECTION, ELSET=EALL, MATERIAL=BRONZE\n\n";
    // The only constraints are transverse displacement constraints on the
    // rotational axis. Axial translation remains free, matching the production
    // axisymmetric model's treatment of the rigid axial mode.
    bool have_axis = false;
    for (const auto& n : mesh.nodes)
        if (std::hypot(n.x,n.y) < 1e-9) { have_axis = true; break; }
    if (have_axis) {
        f << "*NSET, NSET=AXIS\n";
        bool first=true;
        for (const auto& n : mesh.nodes) {
            if (std::hypot(n.x,n.y) < 1e-9) {
                if (!first) f << ", ";
                f << n.id; first=false;
            }
        }
        f << "\n*BOUNDARY\nAXIS, 1, 2, 0.0\n";
    }
    f << "*STEP\n*FREQUENCY, STORAGE=YES\n" << requested_modes << "\n"
      << "*NODE FILE, OUTPUT=3D\nU\n*END STEP\n";
    if (!f) { fail(error, "reference input: write failure"); return false; }
    return true;
}

} // namespace BellOptimizer::Validation
