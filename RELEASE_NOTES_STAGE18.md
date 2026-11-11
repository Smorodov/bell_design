# BellOptimizer 1.2.184 — Stage 18 (2026-11-11)

## Purpose

Address the scale-dependent loss of high-frequency modal content. A fixed modal count gives progressively lower absolute-frequency coverage as a geometrically similar bell grows.

## Implementation

Both `synth-contact` and `synth-clapper` now use adaptive direct-FEM modal evaluation. The initial request is 24 modes and is doubled until the returned spectrum reaches `min(20000 Hz, 0.49 * sample_rate)`, the `--modal-max-modes` cap is reached, or an expanded solve fails. The default cap is 192; accepted values are 24..512. A failed expansion keeps the last successful solution rather than aborting synthesis. Metadata makes incomplete coverage explicit.

## Verification

- GCC 14.2, C++17 Release build: `bell` and three focused test executables built successfully.
- `physical_modal_selection_contract`: passed.
- `clapper_scaling_contract`: passed.
- `clapper_scale_energy_audit`: passed.
- CLI help checked for `synth-contact` and `synth-clapper`.
- `internal_clapper_audio_smoke` was attempted but exceeded the 240-second execution limit; it is not counted as passed.

## Limits

This does not yet prove modal completeness below the ceiling: mesh density, harmonic truncation and eigensolver convergence still constrain the result. Shape/family purity is not implemented. Contact WAV output remains peak-normalized to 0.95, so it does not preserve absolute transferred energy or provide calibrated SPL. The adaptive path still needs a full four-scale audio comparison, mesh/mode convergence study and Windows/MSYS2 validation.
