#include "bell_generator.h"
#include "excitation/meridional_influence.hpp"
#include "ccx/ccx_parser.hpp"
#include "fem/internal_modal_solver.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main() {
    BellMaterial material;
    BellGeometry geometry{658.3175709034034, 355.9758748202059,
                          22.85675347179014, 46.62524667849419,
                          {0.0, 0.130778, 0.216866, 0.45905, 1.0}};
    BellMeshConfig config;
    config.nNodes = 15;
    BellGenerator generator(material, config);
    generator.setGeometry(geometry);
    if (!generator.buildMesh()) return 2;

    BellMaterial soft = material;
    soft.E *= 1e-6;
    BellOptimizer::InternalModalSolver solver(true, 10);
    BellOptimizer::Excitation::MeridionalInfluence modal_influence_1;
    BellOptimizer::ModalProblem problem{generator.getMeshData(), soft, ".", "low_frequency", &modal_influence_1};
    BellOptimizer::ModalData internal;
    if (!solver.solve(problem, internal) || internal.modes.empty()) return 3;
    const double first_hz = internal.modes.front().frequency_hz;
    if (!(first_hz > 0.0 && first_hz < 0.5)) {
        std::cerr << "Expected a valid physical mode below 0.5 Hz, got " << first_hz << " Hz\n";
        return 4;
    }

    const auto path = std::filesystem::temp_directory_path() / "bell_low_frequency_parser_test.dat";
    {
        std::ofstream out(path);
        if (!out) return 5;
        out << " E I G E N V A L U E   O U T P U T\n"
            << " 1 -1.0E-14 0.0 0.0 0.0\n"
            << " 2 2.467401100 1.570796327 0.250000000 0.0\n"
            << " 3 22.20660990 4.712388980 0.750000000 0.0\n"
            << " P A R T I C I P A T I O N   F A C T O R S\n"
            << " 2 1 0 0 0 0 1\n"
            << " 3 1 0 0 0 0 1\n"
            << " E F F E C T I V E   M O D A L   M A S S\n"
            << " 2 1 0 0 0 0 1\n"
            << " 3 1 0 0 0 0 1\n"
            << " T O T A L\n";
    }
    const auto parsed = BellOptimizer::CCXParser::parse_file(path.string(), true);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    if (parsed.modes.size() != 2 || std::abs(parsed.modes[0].frequency_hz - 0.25) > 1e-9) {
        std::cerr << "CCX parser failed to retain a physical sub-0.5-Hz mode; modes="
                  << parsed.modes.size() << '\n';
        return 6;
    }

    std::cout << "low_frequency_mode_policy: PASS\n"
              << "internal_first_mode_hz=" << first_hz << '\n'
              << "ccx_sub_0_5_hz_mode_preserved=1\n";
    return 0;
}
