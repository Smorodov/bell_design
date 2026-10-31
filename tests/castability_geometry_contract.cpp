#include "core/bell_generator.h"
#include <cmath>
#include <iostream>


static BellGeometry make_good() {
    BellGeometry g;
    g.H = 1000.0;
    g.R_mouth = 500.0;
    g.t_top = 10.0;
    g.t_bottom = 30.0;
    g.r_knots = {0.0, 0.15, 0.40, 0.70, 1.0};
    return g;
}

int main() {
    BellMaterial material;
    BellMeshConfig cfg;
    cfg.nNodes = 31;

    {
        BellGenerator gen(material, cfg);
        gen.setGeometry(make_good());
        GeometryValidationResult result;
        if (!gen.buildMesh(&result)) {
            std::cerr << "good geometry rejected: " << result.reason << "\n";
            return 1;
        }
        if (!(result.min_thickness > 0.0 && result.max_thickness >= result.min_thickness)) return 2;
        if (!(result.thickness_ratio < BellConstants::MAX_CAST_THICKNESS_RATIO)) return 3;
    }

    {
        BellGeometry bad = make_good();
        bad.t_top = 2.0;
        bad.t_bottom = 60.0;
        BellGenerator gen(material, cfg);
        gen.setGeometry(bad);
        GeometryValidationResult result;
        if (gen.buildMesh(&result)) {
            std::cerr << "gross thickness variation was accepted\n";
            return 4;
        }
        if (result.reason.find("castability") == std::string::npos) {
            std::cerr << "unexpected rejection reason: " << result.reason << "\n";
            return 5;
        }
    }

    std::cout << "PASS: basic castability rejects only gross wall-thickness pathology.\n";
    return 0;
}
