# BellOptimizer 1.2.184 — Stage 19 (2026-11-11)

## Purpose

Make cross-scale audio comparisons meaningful by allowing contact-driven synthesis to retain its raw numerical amplitude instead of always peak-normalizing the output.

## Changes

- Added `--normalization peak|none` to `synth-contact` and `synth-clapper`.
- `peak` remains the default and preserves legacy peak-normalized PCM16 behavior.
- `none` writes IEEE float32 WAV, avoiding PCM clipping while preserving raw pressure-proxy values.
- Metadata now includes raw acoustic RMS as well as raw peak proxy, normalization gain and output mode.
- Added WAV writer float32 mono/stereo support and contract tests.
- Added reproducible four-scale geometry preparation and WAV/metadata comparison tools.
- Updated the consolidated validation record and scale-audit instructions.

## Verification

GCC 14.2, C++17 Release:
- `bell` built and linked.
- `wav_writer_contract`: passed (PCM16 and float32 mono/stereo).
- `clapper_scale_energy_audit`: passed.
- `physical_modal_selection_contract`: passed.
- `multifidelity_pipeline_contract`: passed.
- CLI help checked for both contact-driven synthesis commands.

## Limitations

- Full four-scale adaptive FEM/audio sweep did not complete within the available execution window; no Stage 19 physical comparison is claimed.
- Float32 raw output preserves the model's relative proxy amplitude but is not calibrated pressure or SPL.
- Full CTest, Windows/MSYS2 and microphone calibration remain outstanding.
