#include "cylindrical_global_assembly.hpp"
#include "cylindrical_q4_element.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BellOptimizer::FEM {
CylindricalGlobalSystem assemble_cylindrical_global_system(
    const BellMeshData& mesh, double E, double nu, double rho,
    const std::vector<int>& harmonics, unsigned angular_samples, double zero_tolerance, double ellipticity) {
    if (!(ellipticity > 0.0 && ellipticity <= 1.0) || !std::isfinite(ellipticity))
        throw std::invalid_argument("cylindrical global assembly: ellipticity must be finite in (0,1]");
    if (mesh.nodes.empty() || mesh.elements.empty())
        throw std::invalid_argument("cylindrical global assembly: mesh is empty");
    if (!(zero_tolerance >= 0.0) || !std::isfinite(zero_tolerance))
        throw std::invalid_argument("cylindrical global assembly: zero tolerance must be finite and nonnegative");
    if (harmonics.empty() || harmonics.front() != 0)
        throw std::invalid_argument("cylindrical global assembly: harmonics must start at zero");
    std::size_t basis_count = 0;
    int previous = -1;
    for (int h : harmonics) {
        if (h < 0 || h <= previous)
            throw std::invalid_argument("cylindrical global assembly: harmonics must be increasing, unique and nonnegative");
        basis_count += (h == 0 ? 1u : 2u);
        previous = h;
    }
    const std::size_t dofs_per_node = 3u * basis_count;
    if (mesh.nodes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / dofs_per_node)
        throw std::length_error("cylindrical global assembly: DOF count exceeds sparse matrix index range");
    const int ndof = static_cast<int>(mesh.nodes.size() * dofs_per_node);
    std::vector<detail::SparseTriplet> kt, mt;
    if (dofs_per_node > std::numeric_limits<std::size_t>::max() / 4u)
        throw std::length_error("cylindrical global assembly: local DOF count overflow");
    const std::size_t local_dofs = 4u * dofs_per_node;
    if (local_dofs != 0 && mesh.elements.size() <= kt.max_size() / local_dofs / local_dofs) {
        const std::size_t reserve = mesh.elements.size() * local_dofs * local_dofs;
        kt.reserve(reserve); mt.reserve(reserve);
    }
    CylindricalGlobalSystem out; out.basis_count = basis_count; out.dofs_per_node = dofs_per_node;
    for (const auto& element : mesh.elements) {
        const int ids[4] = {element.n1, element.n2, element.n3, element.n4};
        std::array<MeridionalNode, 4> local_nodes{};
        for (int a = 0; a < 4; ++a) {
            if (ids[a] <= 0 || static_cast<std::size_t>(ids[a]) > mesh.nodes.size())
                throw std::out_of_range("cylindrical global assembly: element node index out of range");
            local_nodes[static_cast<std::size_t>(a)] = {mesh.nodes[static_cast<std::size_t>(ids[a] - 1)].r,
                                                        mesh.nodes[static_cast<std::size_t>(ids[a] - 1)].z};
        }
        const auto local = assemble_cylindrical_q4(local_nodes, E, nu, rho, harmonics, angular_samples, ellipticity);
        std::vector<int> map(local_dofs);
        for (int a = 0; a < 4; ++a) for (std::size_t c = 0; c < dofs_per_node; ++c)
            map[static_cast<std::size_t>(a) * dofs_per_node + c] =
                static_cast<int>(static_cast<std::size_t>(ids[a] - 1) * dofs_per_node + c);
        for (std::size_t i = 0; i < local_dofs; ++i) for (std::size_t j = 0; j < local_dofs; ++j) {
            const double kval = local.stiffness[i][j], mval = local.mass[i][j];
            if (std::abs(kval) > zero_tolerance) kt.push_back({map[i], map[j], kval});
            if (std::abs(mval) > zero_tolerance) mt.push_back({map[i], map[j], mval});
        }
    }
    out.stiffness = detail::SparseMatrix::from_triplets(ndof, std::move(kt));
    out.mass = detail::SparseMatrix::from_triplets(ndof, std::move(mt));
    // Accumulating adjacent elements in different floating-point orders can
    // leave tiny relative asymmetries in very small off-diagonal entries. The
    // continuum bilinear forms are symmetric, so project the assembled CSR
    // matrices onto their symmetric part before passing them to the eigensolver.
    auto symmetrize = [ndof](detail::SparseMatrix& matrix) {
        std::vector<detail::SparseTriplet> symmetric;
        if (matrix.values.size() <= symmetric.max_size() / 2u)
            symmetric.reserve(matrix.values.size() * 2u);
        for (int row = 0; row < ndof; ++row) {
            for (std::size_t p = matrix.row_offsets[static_cast<std::size_t>(row)];
                 p < matrix.row_offsets[static_cast<std::size_t>(row) + 1u]; ++p) {
                const int col = matrix.columns[p];
                if (col < row) continue;
                const double a = matrix.values[p];
                const double b = matrix.get(col, row);
                const double value = 0.5 * (a + b);
                if (value != 0.0) {
                    symmetric.push_back({row, col, value});
                    if (col != row) symmetric.push_back({col, row, value});
                }
            }
        }
        matrix = detail::SparseMatrix::from_triplets(ndof, std::move(symmetric));
    };
    symmetrize(out.stiffness);
    symmetrize(out.mass);
    return out;
}
} // namespace BellOptimizer::FEM
