#include "commands/optimize.hpp"
#include "commands/synthesize.hpp"
#include "commands/synthesize_geometry.hpp"
#include "commands/synthesize_contact.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#ifndef BELL_VERSION
#define BELL_VERSION "dev"
#endif

#ifdef _WIN32
#include <windows.h>
#endif

namespace {

void print_help() {
    std::cout
        << "BellOptimizer " << BELL_VERSION << "\n\n"
        << "Usage:\n"
        << "  bell optimize [--solver=ccx|internal]\n"
        << "  bell synth [file.dat]\n  bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav]\n  bell synth-contact [geometry.json] [--output file.wav] [--strike-azimuth-deg deg] [--listener-polar-deg deg] [--listener-azimuth-deg deg]\n  bell synth-clapper [geometry.json] [--output file.wav] [--tip-radius m] [--tip-length m] [--initial-angular-velocity rad/s]\n"
        << "  bell --help\n"
        << "  bell --version\n"
        << "\n"
        << "Global option:\n"
        << "  --workdir <dir>    Working directory containing JSON/config/data files\n"
        << "\n"
        << "Examples:\n"
        << "  bell optimize --solver=internal --workdir bells\n"
        << "  bell synth best_bell.dat --workdir bells\n  bell synth-geometry best_geometry.json --ellipticity 0.85 --output bell_elliptic.wav --workdir bells\n  bell synth-contact best_geometry.json --output audio\\bell_contact.wav --workdir .\n  bell synth-clapper best_geometry.json --output audio\\bell_clapper.wav --workdir .\n";
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
        std::cout << "BellOptimizer " << BELL_VERSION << "\n";
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
    if (command == "synth-contact") {
        return command_synthesize_contact(argc - 1, argv + 1);
    }
    if (command == "synth-clapper") {
        return command_synthesize_clapper(argc - 1, argv + 1);
    }
    if (command == "synth" || command == "synthesize") {
        return command_synthesize(argc - 1, argv + 1);
    }

    std::cerr << "Unknown command: " << command << "\n\n";
    print_help();
    return 2;
}
