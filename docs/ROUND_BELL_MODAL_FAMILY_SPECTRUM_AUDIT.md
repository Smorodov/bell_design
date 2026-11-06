# Round-bell modal-family spectrum audit

## Purpose and status

This audit checks whether the direct cylindrical vector-FEM spectrum can safely be mapped to the five acoustic partial families used by the clapper objective. It is a diagnostic, not a calibrated acoustic validation and not a production-objective promotion.

## Reproducible setup

- Geometry: the circular `bells/best_geometry.json` profile (ellipticity fixed at 1).
- Material: the current `BellMaterial` reference preset.
- Meridional mesh: 15 profile nodes, using the production geometry adapter.
- Fourier basis: harmonics `{0, 1, 2, 3}`.
- Angular quadrature: 64 samples.
- Requested modal count: 32; the current evaluator internally requests additional modes for rigid-mode/cutoff safety, so the returned diagnostic spectrum contains 40 modes.
- Minimum frequency: 1 Hz.

The test `clapper_modal_spectrum_family_audit` computes a mass-lumped energy diagnostic for each Fourier basis label in each returned eigenvector. It reports the dominant absolute harmonic and its fraction of the total mass-weighted field energy. It also counts close-frequency pairs using the objective's current relative tolerance (`5e-5`).

## Observations

The first frequency clusters in this coarse diagnostic solve are approximately:

| Cluster | Frequency (Hz) | Modes | Dominant-harmonic mass fraction |
|---:|---:|---:|---:|
| 1 | 123.143 | 2 | 0.587 |
| 2 | 256.252 | 2 | 0.819 |
| 3 | 322.456 | 2 | 0.520 |
| 4 | 442.208 | 2 | 0.827 |
| 5 | 443.457 | 1 | 1.000 |

Across the 40 returned modes, 16 close-frequency candidate pairs had matching dominant absolute harmonic labels, but all 16 were flagged as mixed/ambiguous by the conservative 0.95 dominant-energy threshold. A dominant Fourier label therefore does not provide a reliable physical family identity in this combined solve.

The fifth cluster is also only about 0.28% above the fourth, but is not within the objective's `5e-5` degeneracy tolerance. It must not be merged just because the frequencies are visually close. Conversely, the current frequency-only grouping cannot prove that a cluster contains only one physical family if distinct families happen to coincide numerically.

## Interpretation and limits

1. Frequency rank alone is not a validated mapping to the five intended partials.
2. A dominant Fourier label is not enough to assign family identity in the current combined solve. The measured field-energy fractions show substantial multi-label content; the precise cause (eigenspace mixing, formulation coupling, or both) requires a dedicated operator-level analysis.
3. The 1.2.179 root-sum-square aggregation is invariant to orthogonal basis rotations inside a numerically degenerate group. It does **not** yet establish that the group corresponds to one physical family, nor does it prevent accidental cross-family aggregation.
4. This is a coarse 15-node spectrum. It cannot establish mesh convergence, stable family ordering, or target-partial calibration.

## Required next algorithmic step

Investigate a harmonic-block solve for the circular, axisymmetric model: solve each azimuthal order independently, retain the cosine/sine subspace for each nonzero order, then merge the resulting modal groups by frequency while preserving explicit harmonic-order metadata. This requires carefully generalizing the current assembly and axis-regularity constraints, which presently require the zero harmonic in every solve. Add contracts for the `m=0`, `m=1` rigid/axis constraints and `m>=2` axis constraints before using the block solve in the optimizer.

After harmonic identity is preserved, track radial/axial family identity within each harmonic block, verify near-degenerate subspaces, and only then map families to the five target partials using both frequency and clapper-driven acoustic response. Recheck mesh refinement and candidate-ranking stability before any optimizer integration.

## Production policy

The existing `western_five_partial` objective remains the production default. The clapper acoustic objective remains experimental and must not be used to rank production candidates until harmonic-family identity, final-scale FEM evaluation, clapper parameter calibration, and the documented quantitative gates are satisfied. Chinese/elliptical research remains outside this round-bell task.
