# Consolidated numerical results

This file replaces transient CSV/WAV bundles from research stages 1.2.160–1.2.163. It retains the measured setup, durable numbers and interpretation; generated signals and full sample-by-sample tables are intentionally not shipped.

## Modal/acoustic mesh refinement (1.2.160)

Real FEM, q=0.85, Fourier basis `{0,2,4}`, tracked two-mode family 10:11, receiver distance 3 m, polar angle 0.82 rad, azimuth 0.37 rad, 720 azimuth samples and 16 meridional subdivisions. Profile node counts: 15/23/31/41/49.

| Nodes | Pair frequency change vs 15 nodes | Pair radiation norm change vs 15 nodes |
|---:|---:|---:|
| 15 | 0 | 0 |
| 23 | 2.2256% | 4.0712% |
| 31 | 2.9786% | 5.4739% |
| 41 | 3.3691% | 6.2112% |
| 49 | 3.5233% | 6.5045% |

Adjacent 41→49 changes: frequency 0.1596%, pair radiation norm 0.3127% at each ear. Coarse 15-node result is not converged for this metric. Fixed-profile meridional quadrature 8→16 changed the complex field by 0.241344%. One geometry/family only; not a global acoustic-accuracy claim.

## Finite-distance vs far-field receiver comparison (1.2.161)

Real cylindrical FEM; 13 modal components; both ears; distances 1/2/3/5 m; head yaw 0°/15°. Global absolute complex difference normalized by the strongest modal transfer was about 1.18%. Individual weak-mode relative errors reached about 32% near radiation minima. This is a comparison between approximate kernels, not a physical error estimate.

## Full stereo WAV comparison (1.2.162)

Four matched 10 s, 48 kHz stereo renders used identical FEM/contact/modal release histories at yaw 0° and +15°.

- At 0° yaw, L/R matched to PCM quantization precision.
- At +15°, L/R correlation: 0.9259 far-field; 0.9441 finite-distance.
- Interaural propagation offset: 68.8 μs (~3.30 samples at 48 kHz).
- Finite-distance vs far-field waveform correlation after peak normalization: 0.9489 at 0°; 0.9569/0.9402 for left/right at +15°.
- At 0° yaw, finite-distance vs far-field band-level differences after independent peak normalization: +1.65 dB (50–200 Hz), −0.93 dB (200–500 Hz), +0.61 dB (500–1000 Hz), −0.61 dB (1–2 kHz), −0.88 dB (2–5 kHz), −0.92 dB (5–10 kHz).

No physical reference establishes which path is more accurate. Keep far-field default and finite-distance opt-in.

## Listener geometry matrix (1.2.163)

27 configurations: distance 1/3/5 m, polar angle 0.45/0.82/1.20 rad, yaw 0°/15°/30°, two ears, 11 modes; 594 receiver/mode samples.

- Aggregate normalized complex L2 difference: 0.920958%.
- Maximum eligible per-mode/receiver relative differences: 11.49% at 1 m, 3.78% at 3 m, 2.33% at 5 m.
- Eligibility threshold: far-field amplitude at least 0.1% of the strongest mode in that configuration, to avoid near-zero denominators.

One representative bell and fixed 0.18 m ear spacing. Not a physical error bound.

## Modal tracking closure

- q sequence 0.99→0.96: 4/4 tracker updates accepted; minimum pair similarity 0.999963; no ambiguous updates.
- 15→23 node continuation: update accepted, pair similarity 0.999909; legacy boundary-shape matcher did not report shape tracking.
- Fourier `{0,2,4}`→`{0,2,4,6}`: pair similarity 0.999960.
- Controlled 1.2.159 cross-mesh assignment patch: minimum individual overlap 0.999854; 2D principal similarity 0.999966; 4/4 reference identities assigned one-to-one.

The full-Q4 cross-mesh metric remains separate from the diagnostic Kalman tracker. These results close only the tested diagnostic scope.

## Reproduction policy

The current regression contracts remain in `tests/`. The removed one-shot result generators are not part of the normal build. If a future research question requires regenerating a matrix, add a temporary or separately scoped research executable and document its exact geometry, mesh, modes, quadrature, receiver and normalization before comparing results.
