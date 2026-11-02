# Cylindrical vector-Fourier axis regularity audit — 1.2.134

## Status — implemented in the continuation release

The 1.2.133 fixed-DOF policy was incomplete, and the original audit also had the two `m=1` compatibility signs reversed. Both defects are corrected in the follow-on implementation in `cylindrical_global_modes.*`; the original 1.2.133 archive remains unchanged.

The direct cylindrical Q4 basis stores each displacement component in a real Fourier basis:

- `0`: constant;
- `+m`: `cos(m phi)`;
- `-m`: `sin(m phi)`.

The global DOF layout remains `[node][component][Fourier basis]`, with components `[u_r, u_phi, u_z]`.

## Derivation from Cartesian smoothness

A smooth displacement field at the axis must have finite, direction-independent Cartesian components. The cylindrical-to-Cartesian relations are

\[
 u_x=u_r\cos\phi-u_\phi\sin\phi,\qquad
 u_y=u_r\sin\phi+u_\phi\cos\phi,\qquad
 u_z=u_z.
\]

For the leading transverse field to approach a constant Cartesian vector, its `m=1` traces must have the form

\[
 u_r^{(1)}=a\cos\phi+b\sin\phi,\qquad
 u_\phi^{(1)}=b\cos\phi-a\sin\phi.
\]

Thus, with the stored cosine/sine coefficients,

\[
 u_{r,1c}+u_{\phi,1s}=0,\qquad
 u_{r,1s}-u_{\phi,1c}=0.
\]

The remaining leading trace conditions are:

- `m=0`: `u_r=0` and `u_phi=0` at the axis; `u_z` may be finite.
- `m=1`: the two transverse compatibility relations above; `u_z=0` at the axis.
- `m>=2`: `u_r=u_phi=u_z=0` at the axis.

These are trace conditions at `r=0`. They do not mean that the corresponding fields vanish in a neighborhood of the axis; regular higher-order terms may grow with powers of `r`.

## Defect in the 1.2.132 implementation

The former `production_axis_constraints()` policy fixed only the `m=0` radial trace, left the `m=0` azimuthal trace free, and fixed every nonzero Fourier trace. It therefore removed the admissible `m=1` transverse constant-vector traces rather than enforcing their two compatibility relations. The former fixed-DOF-only solver interface could not express the coupled `m=1` constraints.

## Implemented exact elimination

The evaluator now constructs a constraint transformation `u = P q` at the global-system reduction boundary:

\[
 K_r=P^T K P,\qquad M_r=P^T M P.
\]

At each axis node, `P`:

- omits the constrained `m=0` radial and azimuthal traces;
- retains the `m=0` axial trace;
- when `m=1` is present, retains two independent transverse variables and maps them to the four stored coefficients as `ur(1c)=a`, `uphi(1s)=-a`, `ur(1s)=b`, `uphi(1c)=b`;
- omits the `m=1` axial traces;
- omits all `m>=2` traces.

The implementation forms `P^T K P` and `P^T M P` by sparse triplet accumulation, solves the reduced generalized eigenproblem, and expands eigenvectors with the same map. It does not use penalty stiffness. The builder validates the Fourier layout and rejects a partial `m=1` pair; harmonic families that do not include `m=1` remain valid.

## Regression coverage

`tests/cylindrical_axis_regularity_contract.cpp` verifies the exact reduced mapping, both constant Cartesian transverse translation patterns, elimination of forbidden axis traces, identity mapping away from the axis, and expanded eigenvector compatibility. The circular localized-strike integration test exercises the constrained FEM modes and checks rotation covariance on three profile meshes. A future strengthening item is an independently assembled dense reference for `P^T K P`/`P^T M P` and a full constrained-subspace residual/M-orthogonality audit.

## Scope and non-claims

This derivation concerns smooth vector fields on a circular axis in cylindrical coordinates. It does not validate the current non-circular coordinate mapping at its axis, material behavior, impact-time response, radiation, or stereo. It does not establish frequency convergence. The 1.2.132 three-mesh covariance result remains a provisional historical diagnostic; the follow-on constrained result is recorded below.

## Decision

The exact axis constraint gate is now implemented and its initial contract/integration tests pass. This does **not** authorize production acoustic use: the direct vector-FEM path remains diagnostic, and impact time history plus non-axisymmetric radiation mapping are still missing. Keep the default optimizer and ordinary synthesis on their existing axisymmetric production path.


## Follow-on verification after exact axis elimination

The circular localized-strike integration test now retains the requested number of physical modes after excluding the free-shell rigid-body modes and keeping an additional complete modal pair at the frequency-window boundary. On the tested 15/23/31 profile meshes, the lowest retained pair frequency changed from 81.0672 Hz to 79.2442 Hz (-2.25%), compared with the much larger provisional drift recorded before exact regularity. The lowest-pair generalized-force power changed by -32.6%; this is a material contact/mesh sensitivity and remains diagnostic, not a converged physical result. Maximum rotation-covariance error across the three meshes was 9.11e-5.

The additional mode window is deliberate: six free-body rigid modes can be present when the complete `m=1` basis is included, and a further pair prevents the retained cutoff from splitting a circular cosine/sine family. `requested_modes` continues to mean the number of retained physical modes after the configured minimum-frequency filter.
