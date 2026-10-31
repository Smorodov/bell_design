#pragma once

#include "modal_types.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace BellOptimizer {

struct ModalTrackerConfig {
    // Constant-velocity model in the caller's continuation parameter (for
    // example ellipticity q). Frequencies are measured in Hz.
    double process_frequency_variance = 4.0;
    double process_slope_variance = 25.0;
    double measurement_frequency_variance = 1.0;
    double minimum_subspace_similarity = 0.55;
    double minimum_assignment_margin = 0.035;
    double maximum_log_frequency_error = 0.35;
};

struct ModalTrackEstimate {
    bool valid = false;
    std::size_t candidate_index = 0;
    double predicted_frequency_hz = 0.0;
    double measured_frequency_hz = 0.0;
    double frequency_sigma_hz = 0.0;
    double shape_similarity = 0.0;
    double assignment_score = 0.0;
};

struct ModalStateTrackerReport {
    bool initialized = false;
    bool updated = false;
    bool ambiguous = false;
    bool no_compatible_assignment = false;
    double best_assignment_score = 0.0;
    double pair_subspace_similarity = 0.0;
    double assignment_margin = 0.0;
    std::array<ModalTrackEstimate, 2> tracks{};
};

// A two-track continuation filter. It predicts frequency with a constant-velocity
// Kalman model and associates candidates using sign-invariant full Fourier outer-
// wall traces (or legacy radial/axial traces), resampled to a common normalized
// profile coordinate. Competing assignments are scored using pair-subspace and
// individual-shape continuity. Ambiguous or incompatible assignments are not committed.
class ModalStateTracker {
public:
    explicit ModalStateTracker(const ModalTrackerConfig& config = {});

    bool initialize(const ModalData& modes, std::size_t first_index,
                    std::size_t second_index, double parameter_value);
    ModalStateTrackerReport update(const ModalData& candidates,
                                   double parameter_value);
    void reset();
    bool isInitialized() const { return initialized_; }

private:
    struct Track {
        // State [frequency, d(frequency)/d(parameter)].
        double frequency = 0.0;
        double slope = 0.0;
        double p00 = 1.0, p01 = 0.0, p11 = 1.0;
        std::vector<double> shape;
        std::vector<int> harmonics;
        ModalMode full_mode;
    };

    ModalTrackerConfig config_;
    std::array<Track, 2> tracks_{};
    double last_parameter_ = 0.0;
    bool initialized_ = false;
};

} // namespace BellOptimizer
