#include "fem/cylindrical_global_modes.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace BellOptimizer;
using namespace BellOptimizer::FEM;

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

int dof(int node, int component, int basis, int basis_count) {
    return node * 3 * basis_count + component * basis_count + basis;
}
}

int main() {
    try {
        BellMeshData mesh;
        mesh.nodes.resize(2);
        mesh.nodes[0].r = 0.0;
        mesh.nodes[1].r = 0.1;
        const std::vector<int> harmonics{0, 1, 2};
        constexpr int basis_count = 5; // 0, +1c, -1s, +2c, -2s
        constexpr int dofs_per_node = 3 * basis_count;
        const auto transform = make_cylindrical_axis_regularity_transform(
            mesh, dofs_per_node, harmonics);
        require(transform.reduced_index.size() == 30, "wrong full DOF count");
        require(transform.reduced_dof_count == 18, "wrong reduced DOF count");

        // m=0 transverse traces vanish; axial m=0 survives.
        require(transform.reduced_index[dof(0, 0, 0, basis_count)] == -1,
                "m=0 radial axis trace must be eliminated");
        require(transform.reduced_index[dof(0, 1, 0, basis_count)] == -1,
                "m=0 azimuthal axis trace must be eliminated");
        require(transform.reduced_index[dof(0, 2, 0, basis_count)] >= 0,
                "m=0 axial axis trace must survive");

        // A constant Cartesian x translation has ur=cos(phi), uphi=-sin(phi).
        const int ur1c = dof(0, 0, 1, basis_count);
        const int uphi1s = dof(0, 1, 2, basis_count);
        require(transform.reduced_index[ur1c] == transform.reduced_index[uphi1s],
                "x-translation traces must share one reduced variable");
        require(transform.scale[ur1c] == 1.0 && transform.scale[uphi1s] == -1.0,
                "x-translation signs are incorrect");

        // A constant Cartesian y translation has ur=sin(phi), uphi=cos(phi).
        const int ur1s = dof(0, 0, 2, basis_count);
        const int uphi1c = dof(0, 1, 1, basis_count);
        require(transform.reduced_index[ur1s] == transform.reduced_index[uphi1c],
                "y-translation traces must share one reduced variable");
        require(transform.scale[ur1s] == 1.0 && transform.scale[uphi1c] == 1.0,
                "y-translation signs are incorrect");

        // m=1 axial and all m>=2 traces vanish at the axis.
        for (int b : {1, 2, 3, 4}) {
            require(transform.reduced_index[dof(0, 2, b, basis_count)] == -1,
                    "m>0 axial axis trace must be eliminated");
        }
        for (int c = 0; c < 3; ++c) {
            for (int b : {3, 4})
                require(transform.reduced_index[dof(0, c, b, basis_count)] == -1,
                        "m>=2 axis trace must be eliminated");
        }

        // The non-axis node must remain an identity map.
        for (int c = 0; c < dofs_per_node; ++c) {
            const int full = dof(1, 0, 0, basis_count) + c;
            require(transform.reduced_index[full] >= 0 && transform.scale[full] == 1.0,
                    "non-axis DOF must remain unconstrained");
        }

        // Verify that the exact transformed generalized eigenproblem returns
        // vectors satisfying every axis trace relation after expansion.
        std::vector<detail::SparseTriplet> kt, mt;
        for (int i = 0; i < 30; ++i) {
            kt.push_back({i, i, 1.0 + 0.1 * i});
            mt.push_back({i, i, 1.0});
        }
        CylindricalGlobalSystem system;
        system.stiffness = detail::SparseMatrix::from_triplets(30, std::move(kt));
        system.mass = detail::SparseMatrix::from_triplets(30, std::move(mt));
        system.basis_count = basis_count;
        system.dofs_per_node = dofs_per_node;
        const auto solved = solve_cylindrical_global_modes(
            system, transform, 3, 160, 1e-9, false);
        require(solved.converged && solved.mode_shapes.size() == 3,
                "constrained generalized eigenproblem did not converge");
        for (const auto& mode : solved.mode_shapes) {
            require(std::abs(mode[ur1c] + mode[uphi1s]) < 1e-12,
                    "expanded eigenvector violates m=1 x compatibility");
            require(std::abs(mode[ur1s] - mode[uphi1c]) < 1e-12,
                    "expanded eigenvector violates m=1 y compatibility");
            for (int c = 0; c < 3; ++c) {
                for (int b : {3, 4})
                    require(std::abs(mode[dof(0, c, b, basis_count)]) < 1e-12,
                            "expanded eigenvector has forbidden m>=2 axis trace");
            }
        }
        std::cout << "PASS: cylindrical vector-Fourier axis regularity contract\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
