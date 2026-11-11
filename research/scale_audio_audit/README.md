# Scale-dependent audio level audit

This folder contains reproducible case preparation and analysis tools for comparing round bells at linear scales 0.2, 0.5, 1.0 and 2.0.

## Why this audit exists

The default contact-driven renderers peak-normalize to 0.95. That is useful for auditioning, but removes the absolute amplitude difference between cases. Comparing normalized WAV RMS cannot answer whether the larger bell received less energy or radiated a lower pressure proxy.

Stage 19 adds `--normalization none` to `synth-contact` and `synth-clapper`. It writes IEEE float32 WAV without peak normalization or PCM clipping. The float values preserve the renderer's relative pressure-proxy scale. They are **not calibrated pascals or SPL**, so they support comparisons within the same model/configuration only.

## Prepare four geometries

From the project root:

```sh
python tools/prepare_scale_cases.py --input bells/best_geometry.json --output-dir research/scale_audio_audit/cases
```

Only physical dimensions (`H`, `R_mouth`, `t_bottom`, `t_top`) are scaled. Profile knot fractions and ellipticity remain unchanged.

## Suggested clapper sweep

For each scale `s` in `{0.2, 0.5, 1, 2}`, use:

- mass `5*s^3` kg;
- COM length `0.46*s` m;
- tip length `0.52*s` m;
- tip radius `0.035*s` m;
- initial angular velocity `3.2/sqrt(s)` rad/s;
- contact integration step and impact horizon proportional to `sqrt(s)`.

Example for `s=0.5` (adapt output path as needed):

```sh
bell synth-clapper research/scale_audio_audit/cases/geometry_0.5.json \
  --output research/scale_audio_audit/clapper_0.5_raw.wav \
  --clapper-mass 0.625 --com-length 0.23 --tip-length 0.26 --tip-radius 0.0175 \
  --initial-angular-velocity 4.5254834 --contact-step 1.41421356e-6 \
  --impact-horizon 0.1767767 --normalization none
```

Use identical mesh resolution, sample rate, material preset, listener geometry, damping policy and modal cap for every scale. Preserve each `_metadata.txt` and `_modes.csv` sidecar.

## Compare results

```sh
python tools/compare_scale_audio.py research/scale_audio_audit/clapper_*_raw.wav --csv research/scale_audio_audit/summary.csv
```

The report includes waveform peak/RMS/sum-of-squares, duration, format, energy-transfer metadata, maximum represented modal frequency and per-mode CSV energy sum where available. `--normalization none` should report float32 format (`3/32`).

## Interpretation rules

1. Compare `Initial clapper energy J`, `Bell modal energy at separation J` and transfer fraction before interpreting sound level.
2. Compare raw float32 waveform RMS/peak only when listener distance, radiation model, material, strike location, damping, sample rate and geometry scaling policy are held fixed.
3. Peak-normalized PCM16 outputs are for listening; do not use their RMS as an absolute-energy metric.
4. Modal coverage must be checked for every scale. A cap-hit/incomplete band is not a converged spectrum.
5. Float32 raw proxy output is not calibrated SPL. A microphone calibration and validated radiation/room model are required for physical loudness claims.

The full four-scale adaptive FEM/audio sweep is intentionally not marked complete in Stage 19: an attempted high-resolution FEM run exceeded the execution window. This archive adds level-preserving export and analysis tools so the sweep can be performed reproducibly on the user's Windows/MSYS2 machine or a longer-running build host.
