#include "core/multifidelity_pipeline.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace BellOptimizer::MultiFidelity;

namespace {
void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

bool contains(const std::vector<std::uint64_t>& ids, std::uint64_t id) {
    return std::find(ids.begin(), ids.end(), id) != ids.end();
}

CalibrationEvidence passingCalibration() {
    CalibrationEvidence e;
    e.report_available = true;
    e.representative_domain_covered = true;
    e.full_volume_mass_verified = true;
    e.degenerate_subspaces_verified = true;
    e.objective_semantics_verified = true;
    e.spearman_rank_correlation = 0.95;
    e.top5_recall = 0.98;
    e.top10_recall = 0.97;
    e.maximum_matched_frequency_error = 0.01;
    return e;
}
}

int main() {
    const std::vector<CandidateEvidence> candidates{
        {1, 4.0, MatchStatus::Verified},
        {2, 1.0, MatchStatus::Verified},
        {3, 2.0, MatchStatus::Ambiguous},
        {4, 0.5, MatchStatus::Unmatched},
        {5, 3.0, MatchStatus::FastEvaluationFailed},
        {6, 0.1, MatchStatus::Verified}
    };

    // Default policy must never screen: the entire population reaches stage 2.
    auto decision = selectAccurateStageCandidates(candidates, passingCalibration());
    require(!decision.screening_applied, "screening must be disabled by default");
    require(decision.accurate_stage_candidate_ids.size() == candidates.size(),
            "default path must forward all candidates");

    PromotionPolicy policy;
    policy.enable_screening = true;
    policy.shortlist_size = 1;

    auto weak = passingCalibration();
    weak.top5_recall = 0.80;
    decision = selectAccurateStageCandidates(candidates, weak, policy);
    require(!decision.screening_applied, "failed recall gate must disable screening");
    require(decision.accurate_stage_candidate_ids.size() == candidates.size(),
            "failed calibration must forward all candidates");

    decision = selectAccurateStageCandidates(candidates, passingCalibration(), policy);
    require(decision.screening_applied, "fully accepted calibration should enable screening");
    require(contains(decision.accurate_stage_candidate_ids, 6),
            "best verified fast candidate must be promoted");
    require(contains(decision.accurate_stage_candidate_ids, 3) &&
            contains(decision.accurate_stage_candidate_ids, 4) &&
            contains(decision.accurate_stage_candidate_ids, 5),
            "ambiguous, unmatched and failed candidates must remain fallback");
    require(!contains(decision.accurate_stage_candidate_ids, 1) &&
            !contains(decision.accurate_stage_candidate_ids, 2),
            "only non-shortlisted verified candidates may be screened after promotion");
    require(decision.fallback_candidate_count == 3,
            "fallback count must include every uncertain candidate");

    // Duplicate IDs are emitted only once, including fallback paths.
    const std::vector<CandidateEvidence> duplicates{
        {7, 1.0, MatchStatus::Ambiguous}, {7, 2.0, MatchStatus::Unmatched}
    };
    decision = selectAccurateStageCandidates(duplicates, passingCalibration(), policy);
    require(decision.accurate_stage_candidate_ids.size() == 1,
            "candidate IDs must be unique in the stage handoff");

    std::cout << "multifidelity_pipeline_contract: PASS\n";
    return 0;
}
