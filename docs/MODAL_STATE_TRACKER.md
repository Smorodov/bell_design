# Modal state tracker

## Purpose and status

`core/modal_state_tracker.*` is a diagnostic two-track continuation aid. It is not part of the production optimizer and does not modify the eigensolver's eigenvalues or eigenvectors.

## State model

Each track uses the state

`x = [f, df/ds]ᵀ`,

where `f` is frequency in Hz and `s` is a caller-defined continuation parameter (for example ellipticity). A constant-velocity Kalman model predicts the next frequency and updates its estimate only after an association has been accepted.

## Association

The implementation:

1. extracts outer-wall radial, azimuthal and axial displacement traces, including signed Fourier channels when available;
2. resamples trace features to a common normalized profile coordinate and aligns Fourier channels by signed harmonic label;
3. normalizes features and removes eigenvector-sign ambiguity;
4. evaluates pair-subspace continuity and individual branch continuity;
5. when full modal fields and identical node layouts are available, also evaluates a row-sum lumped-mass weighted full-field similarity;
6. combines shape/subspace evidence with frequency-prediction likelihoods;
7. compares distinct ordered assignments and rejects ambiguous or incompatible cases without committing an update.

## Known limitations

- Boundary traces are registered by normalized profile coordinate, not physical arc length.
- The tracker has not yet been connected to the independent `fem/cylindrical_modal_field.*` cross-mesh comparison utility. Its own full-field feature still uses row-sum lumped mass and matching node layouts.
- The new utility compares complete Q4 fields across meshes using density-weighted integration and a consistent cross-mass product, but currently requires equal ellipticity and density and rejects common-domain coverage below the configured threshold. It is not a persistent state estimator.
- Fourier zero-padding aligns feature dimensions but does not establish convergence of newly introduced harmonic content.
- The tracker retains one committed state per branch. It does not maintain a persistent multi-hypothesis bank or use a particle filter.
- A physically meaningful initial pair must be supplied; the tracker does not discover the intended Chinese-bell pair.
- Large continuation-parameter jumps can violate the local constant-velocity assumption; subdivide the continuation path.

## Validation

`modal_state_tracker_contract` tests frequency-order swaps, sign invariance, boundary sample-count changes, rotations inside a two-dimensional subspace, prediction/update and ambiguity rejection.

`modal_tracker_fem_continuation` exercises real direct-FEM q continuation, profile refinement and Fourier enrichment. The tested q and mesh-refinement paths are accepted; the tested Fourier-enrichment pair is conservatively rejected. This is a safe failure, not a successful enrichment result.

## Next required work

Integrate the new cross-mesh full-field utility into the state tracker only after validating it on real FEM mesh-refinement and Fourier-enrichment sequences. Then consider physical arc-length registration and cross-Gram/principal-angle comparison. Consider persistent multi-hypothesis tracking only if real continuation cases show unacceptable loss from conservative rejection.
