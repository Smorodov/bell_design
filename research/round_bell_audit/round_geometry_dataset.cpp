#include "bell_generator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

struct Case { std::string tag; BellGeometry geometry; };

static BellGeometry makeGeometry(double H, double R, double top, double bottom,
                                 double a, double b, double c) {
    BellGeometry g;
    g.H = H; g.R_mouth = R; g.t_top = top; g.t_bottom = bottom;
    g.r_knots = {0.0, a, b, c, 1.0};
    g.ellipticity = 1.0;
    return g;
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "round_geometry_dataset.csv";
    constexpr std::uint32_t seed = 0xB3112026u;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> hdist(540.0, 780.0);
    std::uniform_real_distribution<double> rdist(300.0, 410.0);
    std::uniform_real_distribution<double> ttopdist(14.0, 29.0);
    std::uniform_real_distribution<double> tbotdist(34.0, 62.0);
    std::uniform_real_distribution<double> adist(0.085, 0.145);
    std::uniform_real_distribution<double> bdist(0.19, 0.29);
    std::uniform_real_distribution<double> cdist(0.43, 0.56);

    std::vector<Case> cases;
    // Explicit boundary-oriented cases precede the seeded interior sample.
    cases.push_back({"edge_compact", makeGeometry(540, 300, 14, 34, 0.085, 0.19, 0.43)});
    cases.push_back({"edge_tall_wide", makeGeometry(780, 410, 29, 62, 0.145, 0.29, 0.56)});
    cases.push_back({"edge_thin_top", makeGeometry(650, 350, 14, 48, 0.11, 0.23, 0.49)});
    cases.push_back({"edge_thick_top", makeGeometry(650, 350, 29, 48, 0.11, 0.23, 0.49)});

    for (int attempt = 0; cases.size() < 64 && attempt < 10000; ++attempt) {
        const double a = adist(rng);
        const double b = std::max(a + 0.065, bdist(rng));
        const double c = std::max(b + 0.12, cdist(rng));
        if (c >= 0.90) continue;
        cases.push_back({"seeded_" + std::to_string(cases.size() - 3),
            makeGeometry(hdist(rng), rdist(rng), ttopdist(rng), tbotdist(rng), a, b, c)});
    }
    if (cases.size() != 64) {
        std::cerr << "Could not create 64 deterministic geometries\n";
        return 1;
    }

    std::ofstream csv(out);
    if (!csv) return 2;
    csv << "index,tag,H_mm,R_mouth_mm,t_top_mm,t_bottom_mm,r1,r2,r3,ellipticity,nodes_41,elements_41,nodes_81,elements_81,min_thickness_mm,max_thickness_mm,status,reason\n";
    csv << std::setprecision(12);
    BellMaterial material;
    int failures = 0;
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& item = cases[i];
        BellMeshConfig cfg41; cfg41.nNodes = 41;
        BellGenerator g41(material, cfg41);
        g41.setGeometry(item.geometry);
        GeometryValidationResult v41;
        bool ok = g41.buildMesh(&v41);
        std::size_t n41 = 0, e41 = 0;
        if (ok) { n41 = g41.getMeshData().nodes.size(); e41 = g41.getMeshData().elements.size(); }
        BellMeshConfig cfg81; cfg81.nNodes = 81;
        BellGenerator g81(material, cfg81);
        g81.setGeometry(item.geometry);
        GeometryValidationResult v81;
        ok = ok && g81.buildMesh(&v81);
        std::size_t n81 = 0, e81 = 0;
        double min_t = 0.0, max_t = 0.0;
        std::string reason;
        if (ok) {
            const auto& m = g81.getMeshData();
            n81 = m.nodes.size(); e81 = m.elements.size();
            min_t = m.min_thickness; max_t = m.max_thickness;
        } else {
            ++failures;
            reason = !v41.passed ? v41.reason : (!v81.passed ? v81.reason : "mesh_generation_failed");
        }
        const auto& g = item.geometry;
        csv << i << ',' << item.tag << ',' << g.H << ',' << g.R_mouth << ',' << g.t_top << ','
            << g.t_bottom << ',' << g.r_knots[1] << ',' << g.r_knots[2] << ',' << g.r_knots[3]
            << ',' << g.ellipticity << ',' << n41 << ',' << e41 << ',' << n81 << ',' << e81
            << ',' << min_t << ',' << max_t << ',' << (ok ? "PASS" : "FAIL") << ',' << reason << '\n';
        std::cout << (ok ? "PASS " : "FAIL ") << i << ' ' << item.tag << '\n';
    }
    const int accepted = static_cast<int>(cases.size()) - failures;
    std::cout << "seed=" << seed << " cases=" << cases.size() << " accepted=" << accepted
              << " mesh_failures=" << failures << " output=" << out << '\n';
    return accepted >= 60 && csv ? 0 : 3;
}
