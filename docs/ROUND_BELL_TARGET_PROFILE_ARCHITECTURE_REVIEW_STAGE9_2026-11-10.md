# Round-bell target-profile architecture and documentation review — Stage 9

Date: 2026-11-10  
Project version: 1.2.184 (unchanged)  
Scope: source/documentation audit and architecture specification only. No production C++ behavior was changed.

## Findings from Stage 8 source

1. `bells/` contains the legacy Western minor-third, Western major-third, neutral and Sysoi-exemplar target files, plus frequency-scaled Western/Orthodox examples at 220/440/880 Hz. There is no arbitrary-length target profile implementation.
2. `TargetData` in `src/core/bell_types.h` uses four-element `std::array` fields for frequency ratios/weights and amplitude ratios/weights.
3. `EvolutionaryOptimizer::computeError` in `src/core/bell_optimizer.cpp` hard-codes frequency slots `{0,2,3,4}` relative to slot `1`, and amplitude slots `{1,2,3,4}` relative to slot `0`.
4. `src/commands/optimize.cpp` copies JSON arrays into those fixed arrays. The CLI checks `tuning_model` against a hard-coded allow-list.
5. The same loader currently reads `tradition` from target JSON to choose the geometry seed provider. This couples acoustic target data to geometry initialization.
6. `src/models/bell_model_profile.cpp` separately contains geometry-family selection and a tuning-model allow-list. Orthodox-specific output labels also remain in the optimizer.
7. `evaluateClapperAcousticObjective` is another experimental objective primitive with fixed four-term mappings. It is not the production candidate objective and must be migrated or explicitly scoped separately when the generic profile API is introduced.

## Documentation inconsistencies corrected

- `docs/TUNING_PROFILES.md` now inventories the exact files present in `bells/`, explains the legacy four-slot limitation and distinguishes acoustic targets from geometry priors.
- `docs/USER_GUIDE.md` and `docs/RUN_OPTIMIZATION_EXAMPLES_RU.md` now explicitly say that arbitrary target counts and arbitrary reference-mode selection are not supported by Stage 8.
- `docs/ARCHITECTURE.md` now records the current implementation caveat and points to the target architecture contract.
- `docs/OPTIMIZER_AUDIT.md` now documents the actual fixed slot mappings and the migration work.
- `docs/DEVELOPMENT_PLAN.md` now orders the generic-profile implementation and regression gates.
- The Stage 8 report has a post-stage correction clarifying that its successful Orthodox run did not demonstrate a generic profile engine.
- The Stage 7 audit is marked historical and superseded where it conflicts with Stage 8.
- `README.md` links to the new architecture specification.

## Intended architecture

`docs/TARGET_PROFILE_ARCHITECTURE.md` defines the selector-based v2 design, explicit frequency/amplitude reference selectors, variable-length target lists, strict validation, v1 compatibility, separation of `seed_family` from acoustic targets, and required acceptance tests. It is a design baseline, not an implemented CLI format.

## Validation performed

- Enumerated the exact `bells/targets*.json` and `bells/evolution*.json` files in the Stage 8 archive.
- Inspected target data structures, JSON loading, production objective, model-profile validation and the separate clapper objective.
- Parsed every `bells/*.json` file with Python's JSON parser.
- Checked Markdown local links and referenced target/evolution file names for existence.
- Checked the packaged archive listing and SHA-256.

No C++ build or CTest run was performed because this stage changes documentation only. The generic target parser/evaluator remains unimplemented.

## Required next implementation stage

Implement generic in-memory target vectors and a strict selector-based v2 loader while preserving v1 conversion. Capture and compare the legacy Western production baseline before changing any objective path. Remove profile-name objective gating and move geometry seed selection to evolution configuration only after regression coverage is in place. Do not change the production default or version until the parity and validation gates pass.
