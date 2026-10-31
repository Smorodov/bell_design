#include "validation/nonaxisymmetric_reference.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace {

BellMeshData make_annular_validation_mesh() {
    BellMeshData mesh;
    mesh.nodes = {
        {1, 100.0, 0.0}, {2, 120.0, 0.0},
        {3, 120.0, 200.0}, {4, 100.0, 200.0}
    };
    mesh.elements = {{1, 2, 3, 4}};
    return mesh;
}

} // namespace

int main() {
    const BellMeshData mesh = make_annular_validation_mesh();
    constexpr double q = 0.75;
    constexpr int nphi = 16;
    const auto ref = BellOptimizer::Validation::build_reference_mesh(mesh, q, nphi);

    if (ref.nodes.size() != mesh.nodes.size() * static_cast<size_t>(nphi)) {
        std::cerr << "node count mismatch\n"; return 2;
    }
    if (ref.elements.size() != mesh.elements.size() * static_cast<size_t>(nphi)) {
        std::cerr << "element count mismatch\n"; return 3;
    }
    if (!(ref.computational_volume > 0.0) || !(ref.physical_volume > 0.0)) {
        std::cerr << "invalid reference volume\n"; return 4;
    }
    const double volume_ratio = ref.physical_volume / ref.computational_volume;
    if (std::abs(volume_ratio - q) > 1e-12) {
        std::cerr << "volume scaling mismatch: " << volume_ratio << '\n'; return 5;
    }

    std::string error;
    const std::string filename = "nonaxisymmetric_reference_contract.inp";
    if (!BellOptimizer::Validation::write_calculix_c3d8_frequency_input(
            ref, BellMaterial{}, filename, 8, &error)) {
        std::cerr << "failed to write reference input: " << error << '\n'; return 6;
    }
    std::ifstream input(filename);
    std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    std::remove(filename.c_str());
    if (text.find("TYPE=C3D8") == std::string::npos ||
        text.find("*FREQUENCY") == std::string::npos ||
        text.find("*NODE") == std::string::npos) {
        std::cerr << "reference input contract mismatch\n"; return 7;
    }

    // The first and last azimuthal sectors must close periodically.
    const auto& e0 = ref.elements.front();
    const auto& elast = ref.elements.back();
    if (e0.n[0] != 2 || e0.n[4] != mesh.nodes.size() + 2 ||
        elast.n[4] != 2 || elast.n[5] != 1) {
        std::cerr << "periodic seam mismatch\n"; return 8;
    }

    // Check the C3D8 center Jacobian sign for the first cell.
    const auto node_at = [&](int id) -> const auto& { return ref.nodes[static_cast<size_t>(id - 1)]; };
    const double dN[8][3] = {
        {-0.125,-0.125,-0.125},{0.125,-0.125,-0.125},{0.125,0.125,-0.125},{-0.125,0.125,-0.125},
        {-0.125,-0.125,0.125},{0.125,-0.125,0.125},{0.125,0.125,0.125},{-0.125,0.125,0.125}};
    double J[3][3]{};
    for (int a=0;a<8;++a) {
        const auto& n=node_at(e0.n[a]);
        J[0][0]+=dN[a][0]*n.x; J[0][1]+=dN[a][0]*n.y; J[0][2]+=dN[a][0]*n.z;
        J[1][0]+=dN[a][1]*n.x; J[1][1]+=dN[a][1]*n.y; J[1][2]+=dN[a][1]*n.z;
        J[2][0]+=dN[a][2]*n.x; J[2][1]+=dN[a][2]*n.y; J[2][2]+=dN[a][2]*n.z;
    }
    const double detJ = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
                      - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
                      + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    if (!(detJ > 0.0) || !std::isfinite(detJ)) {
        std::cerr << "non-positive C3D8 Jacobian: " << detJ << '\n'; return 9;
    }

    std::cout << "nonaxisymmetric_reference_mesh_contract: PASS\n";
    std::cout << "q=" << q << ", azimuthal_divisions=" << nphi
              << ", volume_ratio=" << volume_ratio << '\n';
    return 0;
}
