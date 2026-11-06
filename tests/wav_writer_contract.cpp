#include "audio/wav_writer.hpp"

#include <cmath>
#include <cstdint>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
std::uint16_t u16(const std::vector<unsigned char>& b, std::size_t p) {
    return static_cast<std::uint16_t>(b.at(p) | (static_cast<unsigned>(b.at(p + 1)) << 8u));
}
std::uint32_t u32(const std::vector<unsigned char>& b, std::size_t p) {
    return static_cast<std::uint32_t>(b.at(p)) |
        (static_cast<std::uint32_t>(b.at(p + 1)) << 8u) |
        (static_cast<std::uint32_t>(b.at(p + 2)) << 16u) |
        (static_cast<std::uint32_t>(b.at(p + 3)) << 24u);
}
std::vector<unsigned char> readAll(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}
void checkHeader(const std::vector<unsigned char>& bytes, std::uint16_t channels,
                 std::uint32_t frames, std::uint32_t sampleRate) {
    require(bytes.size() == 44u + frames * channels * 2u, "unexpected WAV file size");
    require(std::string(bytes.begin(), bytes.begin() + 4) == "RIFF", "missing RIFF signature");
    require(u32(bytes, 4) == bytes.size() - 8u, "incorrect RIFF chunk size");
    require(std::string(bytes.begin() + 8, bytes.begin() + 12) == "WAVE", "missing WAVE signature");
    require(std::string(bytes.begin() + 12, bytes.begin() + 16) == "fmt ", "missing fmt chunk");
    require(u32(bytes, 16) == 16u && u16(bytes, 20) == 1u, "not canonical PCM format");
    require(u16(bytes, 22) == channels, "incorrect channel count");
    require(u32(bytes, 24) == sampleRate, "incorrect sample rate");
    require(u32(bytes, 28) == sampleRate * channels * 2u, "incorrect byte rate");
    require(u16(bytes, 32) == channels * 2u && u16(bytes, 34) == 16u, "incorrect block alignment/bit depth");
    require(std::string(bytes.begin() + 36, bytes.begin() + 40) == "data", "missing data chunk");
    require(u32(bytes, 40) == frames * channels * 2u, "incorrect data chunk size");
}
}

int main() {
    try {
        const auto base = std::filesystem::temp_directory_path() / "belloptimizer_wav_writer_contract";
        const auto monoPath = base.string() + "_mono.wav";
        const auto stereoPath = base.string() + "_stereo.wav";
        const auto invalidPath = base.string() + "_invalid.wav";
        const std::vector<double> mono{-1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
        require(BellOptimizer::Audio::writeMonoPcm16(monoPath, mono, 48000), "mono writer failed");
        const auto monoBytes = readAll(monoPath);
        checkHeader(monoBytes, 1, static_cast<std::uint32_t>(mono.size()), 48000);
        require(u16(monoBytes, 44) == 0x8000u, "-1.0 PCM endpoint incorrect");
        require(u16(monoBytes, 46) == 0xc000u, "negative PCM scaling incorrect");
        require(u16(monoBytes, 48) == 0u, "zero PCM sample incorrect");
        require(u16(monoBytes, 52) == 0x7fffu, "positive PCM endpoint incorrect");
        require(u16(monoBytes, 54) == 0x7fffu, "positive clipping incorrect");

        const std::vector<double> left{0.25, -0.25, 0.0};
        const std::vector<double> right{-0.75, 0.75, 0.5};
        require(BellOptimizer::Audio::writeStereoPcm16(stereoPath, left, right, 44100), "stereo writer failed");
        const auto stereoBytes = readAll(stereoPath);
        checkHeader(stereoBytes, 2, 3, 44100);
        require(u16(stereoBytes, 44) == 8192u && u16(stereoBytes, 46) == static_cast<std::uint16_t>(-24576),
                "stereo interleaving/order incorrect");

        require(!BellOptimizer::Audio::writeMonoPcm16(invalidPath,
                {0.0, std::numeric_limits<double>::quiet_NaN()}, 48000),
                "non-finite sample was accepted");
        require(!std::filesystem::exists(invalidPath), "invalid input created a partial output file");
        require(!BellOptimizer::Audio::writeStereoPcm16(invalidPath, {0.0}, {}, 48000),
                "mismatched stereo channels were accepted");
        require(!BellOptimizer::Audio::writeMonoPcm16(invalidPath, {0.0}, 0),
                "invalid sample rate was accepted");
        require(!BellOptimizer::Audio::writeStereoPcm16(invalidPath, {0.0}, {0.0},
                std::numeric_limits<int>::max()), "overflowing stereo byte rate was accepted");
        const auto uniqueSuffix = std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
        const auto missingDirectoryPath = (std::filesystem::temp_directory_path() /
            ("belloptimizer_missing_" + uniqueSuffix) / "output.wav").string();
        require(!BellOptimizer::Audio::writeMonoPcm16(missingDirectoryPath, {0.0}, 48000),
                "unwritable destination was accepted");

        std::filesystem::remove(monoPath);
        std::filesystem::remove(stereoPath);
        std::cout << "WAV writer PCM16 mono/stereo contract passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "WAV writer contract failed: " << error.what() << '\n';
        return 1;
    }
}
