# Chinese modal topology audit — 1.2.125

Status: **RESEARCH / DIAGNOSTIC ONLY**

## Purpose

The previous modal-family audit showed a numerically stable two-dimensional continuation, but the two tracked members had different dominant azimuthal orders. This audit therefore tests whether the reduced FEM surface field itself supports a stable nodal-meridian topology instead of inferring family identity from frequency rank or subspace continuity.

## Method

The canonical reduced FEM case is unchanged:

- reference bell geometry from the Chinese research audit;
- ellipticity continuation reference `q=0.9999` to candidate `q=0.70`;
- 16-mode window;
- Fourier basis `{0,1,2,3,4}`;
- 64 angular quadrature points for assembly;
- 720 angular samples for the diagnostic surface zero-crossing measurement.

For each tracked candidate, the test reports:

1. dominant Fourier order `m`;
2. dominant-harmonic energy purity;
3. radial surface-field zero-crossing counts at several meridional probe locations;
4. strike projections at the long and short axes.

The zero-crossing count is only a numerical topology diagnostic. It is not converted into a historical `(m,n)` label.

## Result

The tracked pair remains the same numerical continuation pair as in 1.2.124:

- subspace similarity: `0.928164`;
- frequencies: approximately `1238.26 Hz` and `1254.26 Hz`;
- dominant `m`: `1` and `4`;
- dominant-harmonic purity: `0.400471` and `0.430326`.

The radial surface zero-crossing counts are not stable as a single `2m` topology across the meridional probe locations:

- first candidate: `8, 2, 8, 7`;
- second candidate: `6, 6, 6, 4`.

This is consistent with a strongly mixed Fourier surface field rather than a clean single-order azimuthal topology. The two tracked modes therefore cannot be assigned a common physical `m`.

The existing strike diagnostic remains finite and asymmetric, but this does not repair the topology mismatch.

## Conclusion

The current reduced FEM does **not** establish the physical Chinese-bell `(m,n)a/(m,n)b` family on the canonical case. The failure is stronger than the previous observation of unresolved `n`: the candidate pair already fails the necessary common-`m` condition, and the surface nodal topology is not stable enough to infer a clean meridional pattern.

Therefore the Chinese two-tone mechanism remains:

`RESEARCH BLOCKED — physical modal-family identification not supported by the current reduced FEM.`

This result is a model-capability finding, not a solver failure. The numerical eigenproblem and continuation remain valid under their existing contracts.

## Architectural decision

Do not:

- retune the modal tracker to force a common family;
- choose a pair by frequency rank;
- tune ellipticity/contact width to obtain a desired topology;
- introduce an arbitrary nodal-count threshold;
- promote this diagnostic into a production Chinese objective.

If Chinese research is reopened, the next model question should be whether a more complete non-axisymmetric/3-D structural representation is required to reproduce the experimentally documented paired-mode topology.
