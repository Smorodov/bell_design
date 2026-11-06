#include "bell_generator.h"
#include "modal_types.hpp"
#include "ccx/ccx_parser.hpp"
#include <iostream>
#include <string>
#include <cassert>
#include <cmath>

int main() {
    BellMaterial material;
    BellMeshConfig config;
    config.nNodes = 15;

    BellGeometry geometry;
    geometry.H = 1000.0;
    geometry.R_mouth = 500.0;
    geometry.t_top = 10.0;
    geometry.t_bottom = 30.0;
    geometry.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};

    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    assert(generator.validateGeometry().passed);
    assert(generator.buildMesh());
    assert(generator.getMeshData().nodes.size() > 0);
    assert(generator.getMeshData().elements.size() > 0);
    assert(generator.getTotalArcLength() > 0.0);
    assert(std::isfinite(generator.getOuterRadius(500.0)));
    assert(std::isfinite(generator.derivative(500.0)));

    // A valid candidate must produce positive Q4 Jacobians at all Gauss points.
    for (const auto& element : generator.getElements()) {
        const auto& nodes = generator.getNodes();
        const Node& n1 = nodes[element.n1 - 1];
        const Node& n2 = nodes[element.n2 - 1];
        const Node& n3 = nodes[element.n3 - 1];
        const Node& n4 = nodes[element.n4 - 1];
        const double g = 1.0 / std::sqrt(3.0);
        for (double xi : {-g, g}) for (double eta : {-g, g}) {
            const double dxi[4] = {-0.25*(1.0-eta), 0.25*(1.0-eta),
                                   0.25*(1.0+eta), -0.25*(1.0+eta)};
            const double deta[4] = {-0.25*(1.0-xi), -0.25*(1.0+xi),
                                    0.25*(1.0+xi), 0.25*(1.0-xi)};
            const Node q[4] = {n1,n2,n3,n4};
            double j00=0,j01=0,j10=0,j11=0;
            for (int a=0;a<4;++a) {
                j00 += dxi[a]*q[a].r; j01 += dxi[a]*q[a].z;
                j10 += deta[a]*q[a].r; j11 += deta[a]*q[a].z;
            }
            assert(j00*j11-j01*j10 > 0.0);
        }
    }

    const auto parsed = BellOptimizer::CCXParser::parse_file(std::string(BELL_SOURCE_DIR) + "/tests/reference/ccx/best_bell.dat", true);
    assert(!parsed.empty());
    assert(parsed.modes.size() >= 5);
    assert(parsed.modes.front().mode_number == 2);
    assert(std::abs(parsed.modes.front().frequency_hz - 460.5574) < 1e-3);

    BellOptimizer::ModalData data;
    data.modes.push_back({});
    data.modes.front().frequency_hz = 440.0;
    assert(data.size() == 1);
    data.clear();
    assert(data.empty());
    return 0;
}
