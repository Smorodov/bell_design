#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace BellOptimizer {
namespace MultiFidelity {

enum class MatchStatus {
    Verified,
    Ambiguous,
    Unmatched,
    FastEvaluationFailed
};

struct CandidateEvidence {
    std::uint64_t candidate_id = 0;
    double fast_score = 0.0; // Lower is better; meaningful only for Verified matches.
    MatchStatus match_status = MatchStatus::Unmatched;
};

struct CalibrationEvidence {
    bool report_available = false;
    bool representative_domain_covered = false;
    bool full_volume_mass_verified = false;
    bool degenerate_subspaces_verified = false;
    bool objective_semantics_verified = false;
    double spearman_rank_correlation = 0.0;
    double top5_recall = 0.0;
    double top10_recall = 0.0;
    double maximum_matched_frequency_error = 1.0;
};

struct PromotionPolicy {
    // Deliberately false by default. A user must explicitly enable screening
    // after the calibration evidence has been reviewed and accepted.
    bool enable_screening = false;
    std::size_t shortlist_size = 0;
    double minimum_spearman = 0.90;
    double minimum_top5_recall = 0.95;
    double minimum_top10_recall = 0.95;
    double maximum_frequency_error = 0.03;
};

struct PromotionDecision {
    std::vector<std::uint64_t> accurate_stage_candidate_ids;
    bool screening_applied = false;
    std::size_t fallback_candidate_count = 0;
    std::string reason;
};

// Produces the candidate IDs to evaluate with the accurate model. Until every
// calibration gate passes, all candidates are forwarded. After promotion,
// ambiguous, unmatched and fast-failed candidates remain in the accurate-stage
// fallback pool and are never discarded by the fast evaluator.
PromotionDecision selectAccurateStageCandidates(
    const std::vector<CandidateEvidence>& candidates,
    const CalibrationEvidence& calibration,
    const PromotionPolicy& policy = {});

} // namespace MultiFidelity
} // namespace BellOptimizer
