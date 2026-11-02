# Acoustic field and stereo architecture audit — 1.2.132

## Scope

This audit continues the 1.2.127 review with one stricter question:

> Is the existing circular-bell production model sufficient to construct a spatial acoustic field that can later be rotated for stereo without recomputing the structural or acoustic solution?

The answer is **no, not yet**. The missing part is more specific than a generic "field object": the production excitation/radiation path currently collapses the azimuthal degree of freedom before synthesis.

The detached Chinese research module remains outside this audit and outside the production build.

## 1. Current production acoustic contract

The production synthesis path reconstructs each audible mode from a scalar signed coefficient

`C_m = -part_rz * Qd * omega`.

`part_rz` is the generalized meridional participation and `Qd` is the surface volume-displacement integral. The synthesis layer retains only magnitude plus a sign-derived phase of `0` or `pi`.

This is a valid contract for the current monophonic scalar synthesis approximation. It is **not** a directional acoustic field.

The current production result therefore has no retained:

- listener azimuth;
- listener elevation;
- complex directional transfer function;
- azimuthal radiation coefficients;
- phase versus observation angle;
- channel-specific observation state.

## 2. Stronger finding: circular symmetry alone does not create stereo information

A circular/axisymmetric structural mode with an axisymmetric excitation produces an axisymmetric acoustic field. In a harmonic representation this is the `m=0` component.

For `m=0`, rotation gives

`exp(-i m theta) = 1`.

Therefore rotating an exactly axisymmetric production field produces **no audible spatial change**.

This is a critical constraint. A future stereo implementation cannot obtain meaningful left/right angular variation merely by wrapping the current scalar modal coefficient in a rotation API.

Meaningful rotation requires at least one of the following to be retained by the production model:

1. a non-axisymmetric excitation field;
2. a non-axisymmetric modal response;
3. a directional radiation model that produces non-zero azimuthal orders.

This does **not** imply reopening the Chinese-bell research. The same architectural requirement applies to a circular bell with a localized physical strike: the strike position must first be represented as an angular excitation, then projected onto the structural modes, then converted to an acoustic directional field.

## 3. What the existing vector FEM data can and cannot provide

`ModalMode` already has optional Fourier boundary samples and complete direct-vector-FEM modal data. Those arrays are useful structural representations.

They do not automatically constitute an acoustic radiation field.

In particular:

- a structural displacement coefficient is not a far-field pressure coefficient;
- a Fourier displacement order is not by itself an acoustic radiation order;
- storing cosine/sine structural components does not define listener distance, propagation phase, or radiation boundary conditions;
- a scalar `Qd` obtained after azimuthal integration has already discarded non-zero angular information.

Consequently the existing Fourier storage should be reused only as an input to a separately defined acoustic radiation layer, not relabelled as an acoustic field.

## 4. Minimal production representation

The smallest representation that is compatible with the requested rotation architecture is a modal directional transfer function.

For each retained structural mode `k`, retain a complex coefficient per azimuthal order `m`:

`H_k,m(r_obs, theta_obs) = A_k,m(r_obs) exp(i m theta_obs)`.

For a fixed observation radius and elevation this can be reduced to

`H_k,m = a_k,m + i b_k,m`.

A rotation by `theta` is then the pure post-processing operation

`H_k,m -> H_k,m exp(-i m theta)`.

No structural solve is required. No acoustic re-solve is required if the directional coefficients have already been computed for the declared observation regime.

The representation should remain independent of stereo channels:

```text
modal response
      |
      v
base directional acoustic field
      |
      +---- R(theta_L) ---> left observation
      |
      +---- R(theta_R) ---> right observation
```

## 5. Observation regime must be fixed before implementation

The phrase "acoustic field" is otherwise ambiguous. Before adding a production class, one observation regime must be selected.

### Option A — far-field directional response

Retain only the angular radiation transfer at a declared distance/asymptotic regime.

Advantages:

- smallest representation;
- rotation is exact in harmonic space;
- no spatial grid is required;
- cheap stereo rendering.

This is the preferred first production target if the application only needs stereo orientation rather than arbitrary near-field acoustics.

### Option B — near-field field

Retain pressure/velocity over a spatial observation surface or volume.

This is substantially larger and requires a propagation/radiation model with an explicit validity domain. It should not be introduced merely to support stereo.

**Decision:** do not implement Option B for the current stereo objective.

## 6. Excitation contract required for a circular bell

The current meridional influence is axisymmetric. A localized strike must therefore enter as a separate angular excitation contract.

Conceptually:

`strike(theta,z) -> generalized modal coefficient q_k -> acoustic directional coefficients H_k,m`.

The strike representation must preserve its angular phase. If the same physical strike is rotated by `theta_s`, its angular coefficients transform as

`S_m(theta_s) = S_m(0) exp(-i m theta_s)`.

This is the property that makes repeated stereo/strike orientations cheap.

The structural solve may be reused when the bell geometry is unchanged; only the modal excitation coefficients need to be rotated/recombined according to the declared linear model.

No production code should assume that an axisymmetric generalized force can reproduce this behavior.

## 7. Reference-frame contract

The 1.2.127 audit found that `part_rz` and the acoustic first moment do not currently expose one shared explicit origin convention.

This remains a gate.

Before a new directional acoustic field is introduced, the production documentation and code must identify:

- structural excitation origin;
- acoustic reference origin;
- listener coordinate origin;
- positive azimuth direction;
- positive rotation direction;
- phase convention.

The current scalar synthesis coefficient must not be silently reinterpreted as a physically calibrated pressure field until this contract is explicit.

## 8. Production design decision

Do **not** add a generic `rotate_stereo()` function to `BellSynthesizer` yet.

That would make the API look complete while the physical field is still absent.

Instead implement the following in order:

1. define a compact `AcousticDirectionalField` contract;
2. define the far-field observation convention;
3. define the non-axisymmetric excitation input needed to populate it;
4. construct the field once from modal response and radiation coefficients;
5. validate harmonic rotation mathematically;
6. only then connect two rotated observations to stereo synthesis.

The field object must be immutable from the point of view of channel rendering: left/right rotations operate on copies or returned values and never mutate the shared base field.

## 9. Mathematical validation gates

The rotation layer must pass these independent contracts before it is connected to audio:

- `R(0)F = F`;
- `R(2*pi)F = F`;
- `R(a)R(b)F = R(a+b)F`;
- for isolated order `m`, phase shift is exactly `-m*theta`;
- real cosine/sine evaluation agrees with complex harmonic rotation;
- rotation preserves coefficient norm/energy;
- rotating one channel does not mutate the other channel or the stored base field;
- changing only stereo angle does not invoke FEM or radiation-field construction.

These gates validate the representation and operation, not the physical radiation model.

## 10. Current status

**Production acoustic field:** PROVISIONAL / not implemented.

**Stereo rotation:** NOT IMPLEMENTED.

**Current scalar modal synthesis:** retained unchanged.

**Chinese research:** frozen and isolated.

**Next implementation target:** impact time history and complex modal response, followed by a validated non-axisymmetric radiation mapping. The direct vector-FEM path remains diagnostic.

No acoustic threshold is promoted to a production acceptance gate by this document.

## 10. 1.2.129 contract implementation

The first minimal representation is now implemented as `AcousticDirectionalField`. It stores a declared observation condition plus complex coefficients for non-negative azimuthal orders. It provides two operations only:

- `evaluate(theta)` evaluates the retained harmonic field;
- `rotated(theta0)` returns the same base field expressed after an azimuthal rotation using `exp(-i*m*theta0)`.

The class is deliberately not populated by the current scalar acoustic solver. This keeps the mathematical rotation mechanism available without pretending that the existing `Qd/part_rz` coefficient contains directional information.

The production gate remains blocked until a radiation layer defines how FEM displacement and excitation produce the complex coefficients.

### Rotation invariants now automated

The contract test verifies:

- `R(0)F = F`;
- `R(2*pi)F = F`;
- `R(theta0)F(theta) = F(theta-theta0)`;
- `m=0` is rotation invariant;
- every coefficient magnitude is preserved by rotation.

These are mathematical representation tests, not acoustic validation.

## 11. 1.2.130 localized-strike audit

A direct audit of `src/excitation/localized_strike.*` establishes that the current localized strike is uniform in azimuth. Its spatial weight depends on meridional `z` and its surface integration uses the full ring factor `2*pi*r*ds`. The API has no strike azimuth and projects only the axisymmetric `outer_surface_radial/axial` boundary arrays.

Therefore this module cannot populate non-zero acoustic azimuthal orders for the circular production path. It must not be reused as if it were a point strike merely because its name says `localized`.

`ModalMode` already contains optional direct-vector-FEM Fourier boundary samples. These are the appropriate structural data to reuse when a physically defined circular strike projection is introduced. They remain structural data and are not themselves acoustic radiation coefficients.

The structural projection step is now implemented in 1.2.131 as `modalImpulseCircular()`. It accepts strike azimuth and physical meridional/circumferential 1-sigma widths, projects the periodized-Gaussian angular load analytically onto the retained cosine/sine Fourier basis, and returns signed real generalized impulse factors. The default circular optimizer is not yet wired to the direct vector-FEM path, and the projection is not a complex modal response or an acoustic radiation model. See `docs/LOCALIZED_STRIKE_AZIMUTHAL_AUDIT.md` for equations and tests.

The next gate is contact/mesh convergence, followed by a declared impact time history and non-axisymmetric acoustic radiation mapping.


## 12. 1.2.131 circular strike projection

`modalImpulseCircular()` now provides the first structural input capable of retaining strike azimuth on a circular bell. It reuses the direct vector-FEM Fourier boundary samples already stored in `ModalMode`; it does not add a new FEM formulation. The projection has an analytic periodized-Gaussian angular integral and no angular quadrature. A unit test verifies the cosine/sine rotation law and `m=0` invariance, and an integration test evaluates an exactly circular (`q=1`) model through the existing direct vector FEM evaluator before applying the projection.

This does not close the directional-acoustics gate. The default circular optimizer still uses the axisymmetric `InternalModalSolver`, and the new projection returns real signed generalized impulse factors only. Impact phase/time history and radiation coefficients are still missing. The structural Fourier field must not be treated as an acoustic field.


## 13. 1.2.132 direct-FEM Fourier layout correction

The 1.2.131 circular-FEM integration test showed that Fourier boundary arrays existed, but a deeper covariance reconstruction revealed a DOF-layout defect. The element/global assembly order is `[node][component][Fourier basis]`, while the evaluator had read selected entries using `[node][basis][component]`. This corrupted exported structural fields and the generalized-influence projection. The full-field diagnostic tracker also used an inconsistent indexing expression when matching different Fourier basis counts.

Both paths now use the assembly layout. The regression reconstructs the complete returned outer-wall modal response and checks that changing strike azimuth rotates the field spatially. For the tested circular `q=1` model and provisional 16-mode window, the relative covariance error was at most `1.44e-4` over profile meshes `nNodes=15,23,31`. This validates the structural rotation law for that tested configuration only. The lowest pair frequency still changes about `14.9%` between the coarsest and finest mesh, so modal frequency convergence is not established. Contact-width convergence, axis regularity derivation, complex impact response and radiation mapping remain open.

## 14. 1.2.134 exact axis regularity follow-on

The direct vector-FEM modal path now uses exact axis constraint elimination. The Cartesian smoothness conditions for the real `m=1` Fourier traces are `ur(1c)+uphi(1s)=0` and `ur(1s)-uphi(1c)=0`; the earlier audit's displayed signs were corrected. The evaluator retains six free-body rigid modes plus one complete pair in its internal solve window, then returns the requested number of physical modes after the minimum-frequency filter.

On the 15/23/31-node circular strike diagnostic, the retained lowest-pair frequency changed by -2.25%, while the rotation-covariance error remained below `9.11e-5`. The lowest-pair generalized-force power still changed by -32.6%, so contact/mesh convergence is not closed. These results validate structural rotation covariance for the tested model/window only.

The next physical task is now impact time history and complex modal response, followed by a separately derived non-axisymmetric radiation mapping. Structural Fourier displacement coefficients must not be relabelled as acoustic radiation coefficients. No change is made to default synthesis or stereo rendering.


## 1.2.135 modal response handoff

A standalone linear impulse-response evaluator now maps signed generalized
modal impulses to damped modal displacement histories. It is deliberately
downstream of structural projection and upstream of acoustic radiation. The
current implementation uses an ideal impulse and a per-mode T60; it does not
yet represent finite hammer contact or create directional pressure. It does
not populate `AcousticDirectionalField` and does not change stereo status.

## 1.2.136 conservative hammer-contact handoff

`src/excitation/hammer_contact.*` now computes a finite-duration unilateral
Hertz contact between a freely moving striker and a supplied set of linear
structural modes. It retains a force time history and modal displacement
histories, with a total-energy audit. The implementation is conservative and
uses fixed-step RK4; it is not a calibrated hammer, does not include contact
losses, and does not generate acoustic radiation coefficients. The test
contract checks non-negative contact force, finite contact duration, striker
impulse balance and energy conservation under time-step refinement.

This new mechanical layer is not connected to the production audio path. The
modal point-shape values and modal normalization must be supplied consistently.
The structural contact/mesh sensitivity identified in 1.2.134 is still open;
this time-integration test does not resolve it. `AcousticDirectionalField`
remains unpopulated by production synthesis and stereo remains unimplemented.


## 1.2.137 contact-driven mono listening path

The new `synth-contact` path uses the axisymmetric modal volume-displacement coefficient already computed by the FEM solver and evaluates a monopole pressure proxy proportional to modal volume acceleration. It is sufficient to produce a repeatable mono listening file, but it does not populate `AcousticDirectionalField`, provide azimuthal directivity, or unblock stereo. The sound file is peak-normalized, not an absolute pressure prediction.


## 1.2.138 directional mono radiation handoff

The opt-in `synth-contact` command now derives a complex far-field transfer coefficient for each retained structural mode by integrating its normal surface displacement over the outer wall with the acoustic phase factor `exp(-i k n·x)`. The structural modal field is solved with real Fourier harmonics m=0..3, and the strike point is evaluated at a configurable azimuth before Hertz contact integration. The listener direction is configurable by polar angle from the bell axis and azimuth.

The complex transfer coefficient is applied to the analytic modal acceleration, preserving the phase shift associated with the listener direction. This is the first production-path use of directional acoustic phase, but it is still a mono point-listener preview. The integration is a first-order far-field approximation with no calibrated radiation constant, measured material loss, room transfer, or BEM acoustic boundary solution. It does not yet populate the reusable `AcousticDirectionalField` data structure and does not implement stereo channel rotation.

The circular rotational covariance check gives 0.9999999985 waveform correlation when strike and listener azimuth are rotated together by 90 degrees. This supports the coordinate convention, not absolute acoustic accuracy. Spectral mesh convergence and contact-force convergence remain separate gates.
