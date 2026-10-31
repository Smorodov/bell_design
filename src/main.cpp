#include "commands/optimize.hpp"
#include "commands/synthesize.hpp"
#include "commands/synthesize_geometry.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#endif

namespace {

void print_help() {
    std::cout
        << "BellOptimizer 1.2.96\n\n"
        << "Usage:\n"
        << "  bell optimize [--solver=ccx|internal]\n"
        << "  bell synth [file.dat]\n  bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav]\n"
        << "  bell --help\n"
        << "  bell --version\n"
        << "\n"
        << "Global option:\n"
        << "  --workdir <dir>    Working directory containing JSON/config/data files\n"
        << "\n"
        << "Examples:\n"
        << "  bell optimize --solver=internal --workdir bells\n"
        << "  bell synth best_bell.dat --workdir bells\n  bell synth-geometry best_geometry.json --ellipticity 0.85 --output bell_elliptic.wav --workdir bells\n";
}

bool set_workdir(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--workdir") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --workdir requires a directory.\n";
                return false;
            }
            try {
                std::filesystem::current_path(argv[++i]);
            } catch (const std::exception& e) {
                std::cerr << "Error: cannot change working directory: " << e.what() << "\n";
                return false;
            }
        } else if (arg.rfind("--workdir=", 0) == 0) {
            try {
                std::filesystem::current_path(arg.substr(10));
            } catch (const std::exception& e) {
                std::cerr << "Error: cannot change working directory: " << e.what() << "\n";
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    if (argc < 2) {
        print_help();
        return 0;
    }

    const std::string command = argv[1];
    if (command == "--help" || command == "-h" || command == "help") {
        print_help();
        return 0;
    }
    if (command == "--version" || command == "-V" || command == "version") {
        std::cout << "BellOptimizer 1.2.96\n";
        return 0;
    }

    if (!set_workdir(argc, argv)) {
        return 2;
    }

    if (command == "optimize") {
        return command_optimize(argc - 1, argv + 1);
    }
    if (command == "synth-geometry" || command == "synthesize-geometry") {
        return command_synthesize_geometry(argc - 1, argv + 1);
    }
    if (command == "synth" || command == "synthesize") {
        return command_synthesize(argc - 1, argv + 1);
    }

    std::cerr << "Unknown command: " << command << "\n\n";
    print_help();
    return 2;
}
