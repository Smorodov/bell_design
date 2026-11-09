# Round-bell optimizer integration hand-off — Stage 4

Date: 2026-11-08  
Baseline: BellOptimizer 1.2.184 / Stage 3 archive  
Status: **integration boundary prepared and exercised; production optimizer not yet modified**

## What was added

A small research-only orchestration boundary was added:

- `research/round_bell_audit/round_bell_candidate_analysis.hpp`
- `research/round_bell_audit/round_bell_candidate_analysis.cpp`

`analyzeRoundBellCandidate(...)` accepts a `BellMeshData`, the `ModalData` produced from that same direct-FEM solve, a `CircularLocalizedStrikeConfig`, and `ModalFamilyAdapterConfig`. It performs the localized circular-strike projection, validates finite/nontrivial projections, invokes the modal-family adapter, and returns both signed projections and grouped diagnostics in one result object. On any failure it clears the result and returns `false`.

The existing `round_optimizer_integration_probe` now uses this exact boundary rather than manually sequencing strike projection and grouping. This makes the probe a test of the hand-off API intended for the future optimizer call site.

## Exact integration boundary in `EvolutionaryOptimizer::evaluateCandidate()`

The current production path has an important gap:

1. For round geometry (`ellipticity == 1.0`), `evaluateCandidate()` currently calls the configured generic `solver_->solve(...)` path.
2. The direct cylindrical vector FEM evaluator is currently called only in the optional elliptical branch (`ellipticity < 1.0 && use_direct_elliptical_fem`).
3. The new adapter requires the complete direct-FEM modal fields (`fem_displacement`, `fem_mass_lumped_weights`, `fem_node_r/z`, `fem_harmonics`) and matching `BellMeshData` for the strike projection. A generic modal result without those fields cannot be used.

Therefore the production integration must **not** simply call `analyzeRoundBellCandidate()` on whatever `solver_` returns. The next production patch must explicitly choose and validate a direct cylindrical FEM solve for the round-bell path, build the matching mesh/model, and pass the exact same mesh and modal result to the hand-off function.

Recommended insertion point: after the selected modal solve has completed and the final modal set has been sorted/tracked, but before analytical geometry/frequency scaling. Keep the analysis as candidate metadata/diagnostics at first. Do not feed family labels into `computeError()` until family identity has been validated over multiple geometries and mesh resolutions.

## Safe integration sequence

1. Add an explicit opt-in round-bell direct-FEM mode (or an equivalent solver policy) rather than silently replacing the current solver path.
2. For each candidate, construct one cylindrical FEM model and evaluate modes from it.
3. Build the circular strike config from an explicit strike-position policy; do not bury a fixed strike location in the objective.
4. Call `analyzeRoundBellCandidate()` with the same model mesh and modal result.
5. Preserve the returned signed projections and family groups as diagnostics. The group projection norm is invariant to orthogonal rotations/sign flips only within an exactly degenerate, mass-orthonormal subspace; near-degenerate and mesh-varying cases still need robustness evidence.
6. Keep the production `western_five_partial` objective and current frequency-ranked slots unchanged during the first integration pass.
7. Compare old and new candidate scores/results in a shadow or audit mode before enabling any objective change.

## Important technical constraints

- The adapter's grouping rule still uses adjacent relative frequency gaps. Chained gaps can create a broad cluster.
- Harmonic energy shares are per-mode and averaged inside each group; this is not a fully subspace-invariant harmonic decomposition for arbitrary mixed subspaces.
- A dominant azimuthal order is a diagnostic classification, not a Hum/Prime/Tertia/Quinta/Nominal identity.
- The `ModalData` passed to analysis must be frequency sorted and its modal fields must refer to the same mesh ordering as the `BellMeshData`.
- The analysis is currently performed on unscaled FEM data. If its representative frequencies are later displayed after the optimizer's analytic scale correction, they must be scaled consistently; do not mix scaled frequencies with unscaled FEM fields silently.
- Strike position and widths are physical/configuration inputs. The current probe values are test settings, not approved production defaults.

## Validation performed

Environment: Linux, GCC 14.2, C++17, Release, `BUILD_SMOKE_TESTS=ON`.

- CMake configure: PASS (`audit_results_stage4_configure.log`).
- Built `round_optimizer_integration_probe`: PASS.
- Built `modal_family_adapter_contract`: PASS.
- Built `circular_localized_strike_projection`: PASS.
- Targeted CTest: **2/2 PASS** (`audit_results_stage4_ctest.log`).
- End-to-end round geometry probe through the new hand-off API: PASS (`audit_results_stage4_probe.log`).
  - modes: 26
  - modal-family groups: 15
  - maximum absolute strike projection: 0.0660095

The source target `bell_core` was built as a dependency. The full application target, full project build, and full CTest suite were not run. Windows/MSYS2 UCRT64 validation was not run.

## Release and scope decision

- Project version remains **1.2.184**.
- No production optimizer source or objective was changed.
- No family-to-western-partial mapping was added.
- Chinese-bell research and SF3 remain out of scope.

This archive is prepared for a controlled, opt-in integration patch. It is not a claim that direct FEM is already active in the production round-bell optimizer.
