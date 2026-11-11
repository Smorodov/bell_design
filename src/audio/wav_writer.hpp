#pragma once

#include <string>
#include <vector>

namespace BellOptimizer::Audio {

// Writes interleaved little-endian PCM16 RIFF/WAVE. Inputs are normalized
// floating-point samples; values outside [-1, 1] are clipped.
bool writeMonoPcm16(const std::string& filename,
                    const std::vector<double>& samples,
                    int sampleRate);

bool writeStereoPcm16(const std::string& filename,
                      const std::vector<double>& left,
                      const std::vector<double>& right,
                      int sampleRate);

// Writes IEEE-float32 RIFF/WAVE without normalization or clipping. This keeps
// relative pressure-proxy amplitudes for comparative analysis; it is not SPL.
bool writeMonoFloat32(const std::string& filename,
                      const std::vector<double>& samples, int sampleRate);
bool writeStereoFloat32(const std::string& filename,
                        const std::vector<double>& left,
                        const std::vector<double>& right, int sampleRate);

} // namespace BellOptimizer::Audio
