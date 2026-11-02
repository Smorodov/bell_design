# Development plan

## Current production direction

The production path remains the Western/circular/elliptic bell solver and synthesis pipeline.

The Chinese two-tone investigation is frozen as a detached research module under `research/chinese_bell/`. Its negative results and reusable numerical techniques are preserved there.

## Next production work

### 1. Audit the circular-bell calculation path

Review the current circular-bell FEM and modal calculations using the lessons that survived the Chinese research:

- exact surface metrics where applicable;
- physically correct local surface normals for generalized forces;
- mass-weighted full-field modal comparison;
- conservative handling of near-degenerate modal subspaces;
- no use of frequency rank as a substitute for physical identity.

Only changes that improve the generic production mechanism should be promoted.

### 2. Define the production acoustic field before stereo

Stereo is a separate production feature and must not depend on the unresolved Chinese modal-family assumptions.

The audit in `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md` found that the current production path stores only scalar acoustic modal coefficients. That representation is insufficient for physical field rotation.

The next step is therefore **not** stereo panning. It is a minimal far-field directional acoustic-field representation with explicit angular phase/directivity, a declared observation regime, and a non-axisymmetric excitation contract. An exactly axisymmetric `m=0` field cannot produce a rotation-dependent stereo result.

The intended architecture is:

- bell geometry and structural modal solve independent of stereo;
- excitation and modal response defined before acoustic rendering;
- acoustic field computed once and retained in a compact representation;
- left/right channel positions obtained by rotating the retained field;
- no FEM or acoustic re-solve when only the stereo angle changes;
- audio-channel mapping outside the structural solver.

The circular localized-strike projection now exists as an explicit API over direct vector-FEM Fourier boundary data. A 1.2.132 audit corrected inconsistent `[node][component][Fourier basis]` indexing in the evaluator and full-field tracker. Actual circular FEM modal-response rotation covariance now passes for the tested 16-mode window. The same mesh sequence shows the lowest pair frequency is not converged (about 14.9% change from nNodes=15 to 31), so this is not yet a production physical gate. Next investigate axis regularity and complete contact/frequency convergence, then define impact time history and non-axisymmetric radiation. Near-field acoustics remains out of scope.

### 3. Implement the minimal directional-field representation

`AcousticDirectionalField` is now present as a representation-only contract. It stores complex azimuthal coefficients and supports exact harmonic rotation without mutating the shared base field. The class is intentionally not populated by the current scalar acoustic path.

The next physical gate is the radiation mapping that produces these coefficients. Stereo rendering remains downstream of that gate.

### 1.2.130–1.2.132: localized strike and modal-layout audit

The legacy `modalImpulse()` remains a meridional/ring-localized diagnostic and is unchanged. Version 1.2.131 adds `modalImpulseCircular()`, which projects a Gaussian surface patch at `(z, phi)` onto the direct vector-FEM cosine/sine boundary modes. Its angular integral is analytic, and its footprint widths are physical lengths in mesh units. Synthetic rotation contracts and an exact `q=1` direct-vector-FEM integration check pass. Version 1.2.132 corrected the evaluator/tracker Fourier DOF layout and added actual modal-response rotation covariance across three meshes. The frequency sequence is still unconverged, and the projection is not connected to the default circular optimizer, which still uses axisymmetric modes.

### 4. Reuse research selectively

The detached Chinese code is a source of candidate numerical techniques, not a production dependency. Each candidate technique must be revalidated against the circular/elliptic production mechanism before adoption.

## Chinese research reopening condition

Do not resume Chinese parameter sweeps merely to obtain a preferred interval. Reopen only if there is a new physical question, new reference evidence, or a more expressive structural model capable of testing the documented `(m,n)a/(m,n)b` mechanism.


### 1.2.134: exact vector-Fourier axis regularity

Implemented `u=Pq` elimination for the cylindrical vector-FEM axis traces, including the coupled `m=1` transverse relations. The sparse solver now reduces both matrices as `P^T K P` and `P^T M P` and expands eigenvectors through the same map. A dedicated contract test and circular strike covariance regression pass. The physical contact/mesh sensitivity remains open; the direct vector-FEM path is still not connected to default synthesis or directional radiation.


## Contact-driven sound generation gate

Version 1.2.138 upgrades the opt-in geometry-to-WAV path to direct cylindrical vector FEM with harmonics m=0..3, azimuth-localized modal contact, and a complex directional far-field surface integral. It remains a normalized mono listening preview. The next gates are modal-frequency and strike-force convergence, calibrated contact parameters, measured acoustic directivity validation, and only then stereo/room rendering. Do not promote it to calibrated physical prediction until those gates pass.


## Internal clapper model gate (design audit)

The current 1.2.138 contact-driven audio path uses a freely translating point striker and samples the outer bell boundary. It is not a model of the suspended internal clapper. The design audit in `docs/CLAPPER_GEOMETRY_AND_KINEMATICS_AUDIT.md` establishes the next sequence: export inner-wall modal boundary fields from the existing full eigenvectors; implement and test closest-point/signed-gap geometry on the true normal-offset inner profile; validate rigid-pendulum kinematics; then couple equal-and-opposite contact forces and audit total energy. Only after those contracts pass should the production contact-audio command be changed. The outer-wall modal data remains reserved for acoustic radiation.

Do not treat the existing Hertz coefficient as calibrated. Derive it from the spherical head, signed local inner-wall curvature and both materials, and report the regime in which Hertz assumptions are valid. Unknown clapper dimensions should remain explicit inputs or sensitivity parameters, not invented defaults.


## Internal clapper implementation foundation (1.2.139)

Completed: export actual inner-wall modal displacement for every Fourier basis component and add a tested spherical-tip/inner-profile closest-point query with cavity-normal and signed-gap output. The contracts compare exported boundary values directly with the full FEM eigenvector and test inside/outside geometry. Next: smooth-profile contact-normal/curvature evaluation, rigid-pendulum kinematics, effective elastic contact coefficient from head and bell material properties, then coupled action/reaction integration and total-energy audit. Do not route production audio through this geometry utility until those later contracts pass.


## Internal clapper physical model (1.2.140)

Completed first end-to-end internal-clapper path: a planar rigid pendulum, true inner-wall modal fields, signed spherical-head gap, equal-and-opposite Hertz force coupling, first-contact/separation detection, conservative energy audit and post-separation directional modal WAV rendering. The command is `bell synth-clapper`; the earlier `synth-contact` remains an external translating-striker diagnostic. The current Hertz coefficient uses a locally-flat-wall approximation and illustrative default mass properties. Remaining gates are measured clapper/pivot data, inner-wall principal-curvature contact correction, mesh and time-step convergence of the coupled impact, contact damping/plasticity evidence, repeated impacts, and measured acoustic validation. Do not add room acoustics or stereo before those physical gates.
