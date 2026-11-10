# Target profile implementation — Stage 10 (2026-11-10)

## Implemented

- Added generic variable-length `FrequencyTarget` and `AmplitudeTarget` records.
- Added a selector-based v2 JSON loader with explicit frequency reference and explicit amplitude reference modes.
- Implemented `sorted_index` selectors (zero-based sorted modal slots); unsupported selector kinds are rejected.
- Added finite/positive ratio validation, finite/non-negative weight validation, duplicate/self-comparison checks, required-field checks, and at-least-one-active-term validation.
- Retained v1 compatibility by converting its four frequency targets `{0,2,3,4}` relative to mode `1` and four amplitude targets `{1,2,3,4}` relative to mode `0` into generic vectors.
- Replaced production optimizer frequency/amplitude scoring with one generic vector evaluator. Profile metadata no longer selects a different objective or special mode labels.
- Frequency scaling now uses the configured frequency-reference mode.
- Geometry seed family is read from `evolution.seed_family`; legacy target `tradition` is a fallback only when no evolution seed family exists. Conflicting values are rejected.
- Added a runnable v2 example and a loader contract test.

## Verification performed

- CMake Release configuration succeeded.
- `bell` CLI target built successfully with GCC 14.2 / C++17.
- `target_profile_loader_contract` built and passed. It checks legacy slot conversion, variable v2 targets, arbitrary reference selection, amplitude reference selection, unsupported selectors, empty objectives, and loads all production `bells/targets*.json` files.
- Short deterministic internal-solver optimization using a v1 Western profile completed successfully.
- Short deterministic internal-solver optimization using the v2 example completed successfully.
- A deliberately conflicting legacy `targets.tradition` vs `evolution.seed_family` combination was rejected as designed.

## Not yet verified / remaining gates

- Full CTest suite was not run to completion in this stage.
- No Windows/MSYS2 build or `.cmd` launcher test was run.
- No frozen pre-Stage-10 Western candidate corpus was available for old-vs-new numerical score and ranking comparison. The v1 conversion preserves the old slot mapping and scoring expression, but full optimization parity is not yet claimed.
- The production default `bells/targets.json` and its historical objective remain unchanged.
- `clapper_acoustic_objective` is a separate experimental primitive and still has its own legacy fixed four-target interface; it is not called by the production evolutionary candidate scorer. It must be migrated or explicitly retired before that diagnostic API is promoted as a general profile consumer.
- Sorted-index selectors do not resolve physical mode identity or near-degenerate circular mode pairs. Tracked-family and degenerate-subspace selectors remain future work.
- Profile `source` metadata is currently informational and is not semantically validated.

## Files

Core: `src/core/bell_types.h`, `src/core/target_profile_loader.hpp`, `src/core/bell_optimizer.cpp`, `src/commands/optimize.cpp`.
Tests: `tests/target_profile_loader_contract.cpp`.
Example: `examples/target_profiles/custom_v2_example.json`.
