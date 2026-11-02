#include "excitation/internal_clapper_geometry.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BellOptimizer;
using namespace BellOptimizer::Excitation;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    try {
        // Vertical inner wall at r=1. The cavity is to the radial-inside (-r) side.
        const std::vector<Point2D> wall{{1.0, 0.0}, {1.0, 2.0}};
        SphericalTipContactGeometry contact;
        require(find_spherical_tip_contact_geometry(wall, 0.7, 1.0, 0.1, contact),
                "valid spherical tip query rejected");
        require(contact.valid && contact.segment_index == 0, "wrong contact segment");
        require(std::abs(contact.segment_parameter - 0.5) < 1e-12, "wrong segment projection");
        require(std::abs(contact.surface_point.r - 1.0) < 1e-12 &&
                std::abs(contact.surface_point.z - 1.0) < 1e-12, "wrong contact point");
        require(std::abs(contact.outward_normal.r + 1.0) < 1e-12 &&
                std::abs(contact.outward_normal.z) < 1e-12, "normal does not point into cavity");
        require(std::abs(contact.centre_distance - 0.3) < 1e-12 &&
                std::abs(contact.gap - 0.2) < 1e-12, "wrong separated gap");

        require(find_spherical_tip_contact_geometry(wall, 0.7, 1.0, 0.4, contact),
                "overlapping tip query rejected");
        require(std::abs(contact.gap + 0.1) < 1e-12, "overlap gap sign is wrong");

        // The same segment viewed from the exterior must not be reported as an
        // internal clapper contact.
        require(!find_spherical_tip_contact_geometry(wall, 1.3, 1.0, 0.1, contact),
                "exterior-side sphere was accepted as internal contact");
        require(!find_spherical_tip_contact_geometry(wall, 0.7, 1.0, 0.0, contact),
                "zero-radius tip was accepted");
        require(!find_spherical_tip_contact_geometry({{1.0, 0.0}}, 0.7, 1.0, 0.1, contact),
                "underspecified profile was accepted");

        // A slanted segment verifies the local normal orientation and projection.
        const std::vector<Point2D> slanted{{1.0, 0.0}, {2.0, 1.0}};
        require(find_spherical_tip_contact_geometry(slanted, 0.5, 1.5, 0.2, contact),
                "valid slanted-wall contact rejected");
        require(std::abs(contact.outward_normal.r + std::sqrt(0.5)) < 1e-12 &&
                std::abs(contact.outward_normal.z - std::sqrt(0.5)) < 1e-12,
                "slanted-wall cavity normal is wrong");

        std::cout << "internal_clapper_geometry_contract: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "internal_clapper_geometry_contract: FAIL: " << e.what() << '\n';
        return 1;
    }
}
