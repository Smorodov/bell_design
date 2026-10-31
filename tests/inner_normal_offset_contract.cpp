#include "bell_generator.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 31;

    BellGeometry geometry;
    geometry.H = 658.3175709034034;
    geometry.R_mouth = 355.9758748202059;
    geometry.t_top = 22.85675347179014;
    geometry.t_bottom = 46.62524667849419;
    geometry.r_knots = {0.0, 0.130778, 0.216866, 0.45905, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    GeometryValidationResult validation;
    if (!generator.buildMesh(&validation)) {
        std::cerr << "mesh build failed: " << validation.reason << '\n';
        return 1;
    }

    const auto& outer = generator.getOuterProfile();
    const auto& inner = generator.getInnerProfile();
    if (outer.size() != inner.size() || outer.size() != static_cast<size_t>(config.nNodes))
        return 2;

    double max_length_error = 0.0;
    double max_orthogonality_error = 0.0;
    for (size_t i = 1; i < outer.size(); ++i) {
        const double dr = inner[i].r - outer[i].r;
        const double dz = inner[i].z - outer[i].z;
        const double thickness = generator.getThicknessAt(outer[i].z);
        const double slope = generator.derivative(outer[i].z);
        max_length_error = std::max(max_length_error,
            std::abs(std::hypot(dr, dz) - thickness) / thickness);
        max_orthogonality_error = std::max(max_orthogonality_error,
            std::abs(slope * dr + dz) / (thickness * std::hypot(1.0, slope)));
        if (!(inner[i].r > 0.0) || !(inner[i].z > inner[i - 1].z)) return 3;
    }

    std::cout << "max_normal_offset_length_relative_error=" << max_length_error << '\n';
    std::cout << "max_normal_orthogonality_error=" << max_orthogonality_error << '\n';
    return (max_length_error < 1e-10 && max_orthogonality_error < 1e-10) ? 0 : 4;
}
