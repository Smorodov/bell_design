# Validation status and measured results

## Stage 16 — clapper scale and energy audit (2026-11-11)

A same-profile round bell was scaled to heights 0.2, 0.5, 1.0 and 2.0 m. Each run used a 21-node meridional mesh, harmonics `{0,1,2,3}`, the same reference material, and geometrically scaled clapper inputs (`m~s^3`, `I~s^5`, lengths/tip radius `~s`, initial angular velocity `~s^-1/2`). The run used 24 requested modes (the solver retained 32 with its modal-pair buffer), 16 kHz output and 1 s audio duration. Energy drift stayed below 3.4e-6.

| Bell height | Initial clapper energy | Bell modal energy at separation | Transfer fraction | Peak contact force | Highest represented mode | Raw acoustic peak proxy | Normalized WAV RMS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 m | 0.008667 J | 0.001153 J | 13.31% | 433 N | 4261 Hz | 4.94 | 0.2176 |
| 0.5 m | 0.33856 J | 0.04775 J | 14.10% | 4.68 kN | 1704 Hz | 20.02 | 0.2196 |
| 1.0 m | 5.41696 J | 0.79292 J | 14.64% | 28.39 kN | 852 Hz | 58.45 | 0.2175 |
| 2.0 m | 86.6714 J | 13.1019 J | 15.12% | 172.08 kN | 426 Hz | 168.03 | 0.2179 |

Interpretation:

- Under the selected Froude-like scaling, initial energy rises as `s^4`; absolute energy transferred to the bell also rises strongly. The transfer fraction increased modestly from 13.3% to 15.1%, rather than decreasing for the larger bell.
- Raw acoustic peak proxy increased with size in this test. It is not calibrated SPL, but it does not support the hypothesis that the larger bell received less energy in this model.
- Normalized WAV RMS is almost identical because every output is peak-normalized. These files cannot be used to compare absolute loudness.
- The represented upper frequency falls approximately inversely with bell size because `synth-clapper` uses a fixed mode count. At 1 m the highest represented mode was only 852 Hz; at 2 m it was 426 Hz. This is a major unresolved timbre-coverage defect and can plausibly make large bells sound duller or perceptually quieter.
- The simulated first-contact point was at 33.9% of bell height for all four cases, outside the usual 50–70% band. The trajectory/pivot setup needs review before treating these values as a representative real bell strike.

The same audit also added a synthetic similarity contract with scales 0.2, 0.5, 1 and 2 m. It verifies `s^4` initial energy scaling, energy closure and an approximately scale-stable transfer fraction for the synthetic modal basis. Its measured transfer fractions were 2.176e-5, 2.245e-5, 2.290e-5 and 2.329e-5; these are deliberately synthetic contract values and must not be confused with the full FEM sweep above. This is a contract test, not a substitute for the full FEM scale sweep.

## Current focused software checks

GCC 14.2 / C++17 Release build:

- `bell` — built and linked;
- `clapper_dynamics_contract` — passed;
- `clapper_scaling_contract` — passed;
- `clapper_scale_energy_audit` — passed;
- `multifidelity_pipeline_contract` — passed;
- `physical_modal_selection_contract` — passed;

A clean GCC 14.2 / C++17 Release build of `bell` and the standalone synthesis SDK succeeded. Eight focused contracts were built and run together: `synthesis_sdk_link_contract`, `multifidelity_pipeline_contract`, `physical_modal_selection_contract`, `optimizer_shape_scale_audit`, `hammer_contact_contract`, `clapper_dynamics_contract`, `clapper_scaling_contract` and `clapper_scale_energy_audit`: **8/8 passed**. A full all-target build and full CTest suite have not been run for this checkpoint.

## Existing numerical evidence retained

- Stage 13 full consistent-mass `m=0` cross-match: 305/305 matches across 61 geometries; minimum overlap 0.999999999999; maximum relative frequency difference 8.29e-13.
- Stage 13 near-degenerate rank-2 subspace audit: minimum principal overlap 0.999999963871. Individual vector correspondence can be lower because bases rotate inside a near-degenerate subspace.
- Stereo/radiation and modal tracking evidence is summarized in `research/CONSOLIDATED_NUMERICAL_RESULTS.md`.

These are scoped numerical results, not blanket physical validation.

## Open gates

1. Increase `synth-clapper` modal coverage by absolute-frequency target, not a fixed count; report incomplete coverage explicitly and control solver cost.
2. Evaluate modal-family/shape purity with a calibrated subspace metric; do not infer it from acoustic amplitude or eigensolver residuals.
3. Revisit the clapper trajectory because the present scale sweep strikes at 33.9% height.
4. Calibrate clapper mass, pivot, tip and initial motion against measured hardware before interpreting transfer fractions as real-world values.
5. Validate mesh, mode-count and time-step convergence of transfer fraction and acoustic output.
6. Establish multi-fidelity shortlist recall before enabling candidate rejection by default.
7. Run full all-target build/CTest and Windows/MSYS2 validation before release promotion.


## Stage 18 — adaptive absolute-frequency modal coverage

The `synth-contact` and `synth-clapper` command paths now retry the direct cylindrical FEM solve with an increasing requested mode count (initial request 24, doubling to the configured `--modal-max-modes` cap; default 192). Expansion stops when the returned spectrum reaches `min(20000 Hz, 0.49 * sample rate)`, the cap is reached, or a larger solve fails. A failed expansion retains the last successful result and is reported in metadata. This is a practical coverage mechanism, not a guarantee that every mode below the ceiling has been found: completeness still depends on mesh resolution, solver convergence and harmonic truncation.

Focused Stage 18 verification: GCC 14.2 Release C++17 build of `bell`, `physical_modal_selection_contract`, `clapper_scaling_contract`, and `clapper_scale_energy_audit`; the three focused tests passed. CLI help was checked for both synthesis commands. A full CTest run, full scale sweep using the adaptive path, Windows/MSYS2 validation, measured SPL calibration and physical modal-family purity validation remain outstanding.

## Stage 19 — preserve raw synthesis level for scale comparisons

Both `synth-contact` and `synth-clapper` accept `--normalization peak|none` (default `peak`, preserving prior behavior). `peak` writes peak-normalized PCM16 as before. `none` writes IEEE float32 WAV without gain normalization or clipping, retaining the synthesis pressure-proxy scale. Metadata records raw acoustic peak/RMS proxies, the applied gain and the output format policy. These values are not calibrated SPL.

Added `tools/prepare_scale_cases.py` to generate geometrically similar round-bell fixtures at 0.2/0.5/1/2 linear scale, and `tools/compare_scale_audio.py` to summarize WAV peak/RMS/sum-of-squares alongside energy-transfer, modal-band and per-mode CSV metadata. Detailed usage and comparison constraints are in `research/scale_audio_audit/README.md`.

Verification: GCC 14.2 / C++17 Release `bell` build passed; `wav_writer_contract`, `clapper_scale_energy_audit`, `physical_modal_selection_contract` and `multifidelity_pipeline_contract` passed (4/4). Both CLI help pages expose `--normalization`. The full four-scale adaptive FEM audio sweep was attempted but exceeded the available execution window; no new physical scale-comparison result is claimed. Windows/MSYS2 and full CTest remain open.
