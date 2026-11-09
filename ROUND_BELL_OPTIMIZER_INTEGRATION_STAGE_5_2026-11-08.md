# Round-bell optimizer integration — Stage 5 (2026-11-08)

## Purpose

Stage 5 connects the direct cylindrical FEM + circular localized-strike + modal-family adapter to the real `EvolutionaryOptimizer::evaluateCandidate()` execution path as an **opt-in diagnostic side channel**.

It intentionally does not replace the current production modal solver, does not change the `western_five_partial` objective, does not change modal tracking, and does not assign FEM modes to Hum/Prime/Tertia/Quinta/Nominal.

## Implementation

- Added `EvolutionParams::enable_round_fem_diagnostics` (default `false`).
- Added `EvolutionParams::round_fem_diagnostics_csv` (default `round_fem_diagnostics.csv`).
- Both options are read from `evolution.json` by `src/commands/optimize.cpp`.
- When the option is enabled and `ellipticity == 1` (within `1e-12`), `evaluateCandidate()` runs a direct cylindrical FEM solve with harmonics `{0,1,2,3}`, then invokes `Research::analyzeRoundBellCandidate()` using the circular strike configuration.
- A compact CSV row is appended with geometry, number of FEM modes, family count, and maximum absolute strike projection.
- Any diagnostic failure is reported to stderr and does not invalidate or modify the production candidate.
- The diagnostic code runs before the existing scale-to-target and objective calculation. Its results are never assigned to `raw_data` or `out_data`.
- Adapter implementation files are linked into `bell_core` so the optimizer can call the analysis boundary.

## Enabling it

Add these keys to the active `bells/evolution.json` (or the selected work directory's `evolution.json`):

```json
{
  "enable_round_fem_diagnostics": true,
  "round_fem_diagnostics_csv": "round_fem_diagnostics.csv"
}
```

The default is off. Remove the keys or set the flag to `false` to restore the previous execution path without diagnostic FEM overhead.

## Validation performed

Environment: Linux, GCC 14.2, C++17, Release.

- CMake configure: PASS.
- `bell_core` build: PASS.
- `bell` CLI build: PASS.
- `modal_family_adapter_contract`: PASS.
- `circular_localized_strike_projection`: PASS.
- Targeted CTest: 2/2 PASS.
- Standalone `round_optimizer_integration_probe`: PASS; 26 modes and 15 groups for its fixture.
- End-to-end CLI smoke run: `bell optimize --solver=internal` with population 2, one generation, deterministic seed 17, and diagnostic flag enabled: completed successfully.
- Runtime diagnostic CSV contained two candidate records:
  - 26 modes, 16 groups, max absolute projection `0.114671139057`.
  - 26 modes, 15 groups, max absolute projection `0.0787861844317`.
- The runtime test retained the legacy optimizer's objective/solver path; diagnostic FEM results were written separately.

Logs and the exact runtime test configuration are included in this archive.

## Not performed / remaining limitations

- Full CTest suite was not run.
- Windows/MSYS2 UCRT64 build was not run.
- Diagnostic side channel is opt-in and adds a second FEM solve per eligible candidate; it is not intended to be enabled for large production searches without measuring runtime.
- Current optimizer diagnostic CSV is a compact per-candidate summary. Per-mode and per-family detail remains available through the research probe CSV outputs.
- FEM modal families are not yet equivalent to the five target partials. No target-slot mapping is approved.
- No production solver replacement or fitness contribution is made in this stage.

## Gate for Stage 6

Before enabling direct FEM as the production evaluator, perform paired-candidate comparisons on a controlled round-bell geometry set: legacy solver vs direct FEM, mesh refinement, frequency ordering and mode-family stability, and a separately validated mapping to the five partial target model. Keep this diagnostic flag disabled by default until those comparisons pass.

Version remains **1.2.184**. This is an integration-stage archive, not a functional release.
