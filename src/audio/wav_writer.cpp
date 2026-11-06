#include "wav_writer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>

namespace BellOptimizer::Audio {
namespace {

void writeU16(std::ofstream& file, std::uint16_t value) {
    const char bytes[2] = {
        static_cast<char>(value & 0xffu),
        static_cast<char>((value >> 8u) & 0xffu)
    };
    file.write(bytes, 2);
}

void writeU32(std::ofstream& file, std::uint32_t value) {
    const char bytes[4] = {
        static_cast<char>(value & 0xffu),
        static_cast<char>((value >> 8u) & 0xffu),
        static_cast<char>((value >> 16u) & 0xffu),
        static_cast<char>((value >> 24u) & 0xffu)
    };
    file.write(bytes, 4);
}

bool validSamples(const std::vector<double>& samples) {
    return !samples.empty() &&
        std::all_of(samples.begin(), samples.end(),
                    [](double value) { return std::isfinite(value); });
}

std::int16_t toPcm16(double sample) {
    const double clipped = std::clamp(sample, -1.0, 1.0);
    const double scaled = clipped < 0.0 ? clipped * 32768.0 : clipped * 32767.0;
    return static_cast<std::int16_t>(std::lrint(scaled));
}

bool writeHeader(std::ofstream& file, std::uint32_t dataBytes,
                 std::uint16_t channels, std::uint32_t sampleRate) {
    const std::uint32_t blockAlign = static_cast<std::uint32_t>(channels) * 2u;
    const std::uint64_t byteRate64 = static_cast<std::uint64_t>(sampleRate) * blockAlign;
    if (byteRate64 > std::numeric_limits<std::uint32_t>::max()) return false;

    file.write("RIFF", 4);
    writeU32(file, 36u + dataBytes);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    writeU32(file, 16u);
    writeU16(file, 1u); // PCM
    writeU16(file, channels);
    writeU32(file, sampleRate);
    writeU32(file, static_cast<std::uint32_t>(byteRate64));
    writeU16(file, static_cast<std::uint16_t>(blockAlign));
    writeU16(file, 16u);
    file.write("data", 4);
    writeU32(file, dataBytes);
    return static_cast<bool>(file);
}

bool validSampleRate(int sampleRate, std::uint32_t channels) {
    if (sampleRate <= 0) return false;
    const std::uint64_t bytesPerFrame = static_cast<std::uint64_t>(channels) * 2u;
    return static_cast<std::uint64_t>(sampleRate) <=
        std::numeric_limits<std::uint32_t>::max() / bytesPerFrame;
}

bool validDataSize(std::size_t frames, std::uint32_t channels) {
    constexpr std::uint64_t maxData =
        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) - 36u;
    const std::uint64_t bytesPerFrame = static_cast<std::uint64_t>(channels) * 2u;
    return static_cast<std::uint64_t>(frames) <= maxData / bytesPerFrame;
}

bool finishFile(std::ofstream& file) {
    file.flush();
    if (!file) return false;
    file.close();
    return !file.fail();
}

} // namespace

bool writeMonoPcm16(const std::string& filename,
                    const std::vector<double>& samples, int sampleRate) {
    if (!validSampleRate(sampleRate, 1u) || !validSamples(samples) || !validDataSize(samples.size(), 1u))
        return false;
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file) return false;
    const auto dataBytes = static_cast<std::uint32_t>(samples.size() * 2u);
    if (!writeHeader(file, dataBytes, 1u, static_cast<std::uint32_t>(sampleRate))) return false;
    for (double sample : samples) writeU16(file, static_cast<std::uint16_t>(toPcm16(sample)));
    return finishFile(file);
}

bool writeStereoPcm16(const std::string& filename,
                      const std::vector<double>& left,
                      const std::vector<double>& right, int sampleRate) {
    if (!validSampleRate(sampleRate, 2u) || left.size() != right.size() ||
        !validSamples(left) || !validSamples(right) || !validDataSize(left.size(), 2u))
        return false;
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file) return false;
    const auto dataBytes = static_cast<std::uint32_t>(left.size() * 4u);
    if (!writeHeader(file, dataBytes, 2u, static_cast<std::uint32_t>(sampleRate))) return false;
    for (std::size_t i = 0; i < left.size(); ++i) {
        writeU16(file, static_cast<std::uint16_t>(toPcm16(left[i])));
        writeU16(file, static_cast<std::uint16_t>(toPcm16(right[i])));
    }
    return finishFile(file);
}

} // namespace BellOptimizer::Audio
