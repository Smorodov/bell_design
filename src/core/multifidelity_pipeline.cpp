#include "multifidelity_pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <utility>

namespace BellOptimizer {
namespace MultiFidelity {
namespace {

bool calibrationPasses(const CalibrationEvidence& e, const PromotionPolicy& p) {
    return e.report_available &&
           e.representative_domain_covered &&
           e.full_volume_mass_verified &&
           e.degenerate_subspaces_verified &&
           e.objective_semantics_verified &&
           std::isfinite(e.spearman_rank_correlation) &&
           e.spearman_rank_correlation >= p.minimum_spearman &&
           std::isfinite(e.top5_recall) && e.top5_recall >= p.minimum_top5_recall &&
           std::isfinite(e.top10_recall) && e.top10_recall >= p.minimum_top10_recall &&
           std::isfinite(e.maximum_matched_frequency_error) &&
           e.maximum_matched_frequency_error <= p.maximum_frequency_error;
}

void appendUnique(PromotionDecision& decision,
                  std::unordered_set<std::uint64_t>& seen,
                  std::uint64_t id) {
    if (seen.insert(id).second) decision.accurate_stage_candidate_ids.push_back(id);
}

} // namespace

PromotionDecision selectAccurateStageCandidates(
    const std::vector<CandidateEvidence>& candidates,
    const CalibrationEvidence& calibration,
    const PromotionPolicy& policy) {
    PromotionDecision decision;
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(candidates.size());

    if (candidates.empty()) {
        decision.reason = "No candidates were supplied.";
        return decision;
    }

    if (!policy.enable_screening || policy.shortlist_size == 0 ||
        !calibrationPasses(calibration, policy)) {
        for (const auto& candidate : candidates)
            appendUnique(decision, seen, candidate.candidate_id);
        decision.reason = !policy.enable_screening
            ? "Screening is disabled; all candidates are forwarded."
            : (policy.shortlist_size == 0
                ? "Shortlist size is zero; all candidates are forwarded."
                : "Calibration gates are not all satisfied; all candidates are forwarded.");
        return decision;
    }

    std::vector<const CandidateEvidence*> ranked;
    ranked.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.match_status == MatchStatus::Verified &&
            std::isfinite(candidate.fast_score)) {
            ranked.push_back(&candidate);
        } else {
            appendUnique(decision, seen, candidate.candidate_id);
            ++decision.fallback_candidate_count;
        }
    }

    std::stable_sort(ranked.begin(), ranked.end(), [](const auto* a, const auto* b) {
        if (a->fast_score != b->fast_score) return a->fast_score < b->fast_score;
        return a->candidate_id < b->candidate_id;
    });

    const std::size_t keep = std::min(policy.shortlist_size, ranked.size());
    for (std::size_t i = 0; i < keep; ++i)
        appendUnique(decision, seen, ranked[i]->candidate_id);

    decision.screening_applied = true;
    decision.reason = "Calibration gates passed; verified candidates were shortlisted and uncertain candidates retained as fallback.";
    return decision;
}

} // namespace MultiFidelity
} // namespace BellOptimizer
