# Localized strike azimuthal audit — 1.2.132

## Conclusion

The current `src/excitation/localized_strike.*` implementation is **meridional/ring-localized**, not a point strike on the circular wall.

Its Gaussian weighting is a function of the meridional coordinate `z` only and the surface integration is explicitly multiplied by `2*pi*r`. Therefore the load is uniform in azimuth:

`j(z,phi) = j(z)`.

Its modal projection is consequently an axisymmetric projection. The function receives only `ModalMode::outer_surface_radial` and `outer_surface_axial`, which are the axisymmetric boundary views. The optional Fourier boundary arrays are not used.

This is not a bug in the current diagnostic itself; it is a hard architectural limitation for directional acoustics.

## Consequence for stereo

A circular bell can have non-axisymmetric structural modes and a physically localized strike can excite them. The present `modalImpulse()` API cannot represent that excitation because it has no strike azimuth and no angular load basis.

Therefore it cannot be used as the production source of non-zero acoustic azimuthal orders `m > 0`.

The correct future chain is:

```text
physical strike position (z, phi)
        |
        v
surface load j(z, phi)
        |
        v
Fourier load coefficients S_m(z)
        |
        v
modal generalized coordinates q_k
        |
        v
radiation mapping
        |
        v
AcousticDirectionalField
```

For a rotation of the same physical strike by `theta_s`, the angular load coefficients must transform consistently with the declared convention. With the field convention used by `AcousticDirectionalField`, the coefficient transformation is:

`S_m(theta_s) = S_m(0) * exp(-i*m*theta_s)`.

## Existing data that can be reused

`ModalMode` already retains, for direct vector FEM, optional outer-surface Fourier components:

- `outer_surface_harmonics`;
- `outer_surface_fourier_radial`;
- `outer_surface_fourier_azimuthal`;
- `outer_surface_fourier_axial`.

These are **structural modal boundary data**, not acoustic radiation coefficients. They are sufficient to investigate a modal projection of an angular strike without adding another FEM solve, provided the excitation is defined consistently with the same Fourier basis.

## Important distinction

The current production circular acoustic coefficient

`C_k = -part_rz * Qd * omega`

comes from an axisymmetric generalized rotation and an azimuthally integrated surface displacement `Qd`. It cannot be retroactively decomposed into directional coefficients.

Likewise, `MeridionalInfluence` is intentionally an axisymmetric generalized influence and must not be relabelled as a localized point-strike model.

## 1.2.131: circular localized-strike projection contract

A separate overload is now implemented as `modalImpulseCircular()` in
`src/excitation/localized_strike.*`. It does not replace or change the legacy
ring-localized `modalImpulse()`.

The new configuration declares:

- `center_z`: axial center of the meridional Gaussian;
- `meridional_sigma`: 1-sigma axial width;
- `circumferential_sigma`: 1-sigma physical arc-length width;
- `azimuth_rad`: center angle of the strike;
- normalized normal/tangent direction fractions.

Both widths use the same length unit as the supplied mesh. In the direct vector
FEM path this mesh is SI, so the caller supplies metres. The circumferential
width is converted locally to angular width `sigma_phi = sigma_s / r` at each
surface segment. It is a declared contact-footprint parameter, not a measured
mallet dimension or an optimizer variable.

The periodized Gaussian angular load is projected analytically. For the
Fourier basis labels `0`, `+m` (cosine) and `-m` (sine), its normalized angular
moments are:

`W_0 = 1`

`W_+m = exp(-m^2 sigma_phi^2 / 2) cos(m phi_s)`

`W_-m = exp(-m^2 sigma_phi^2 / 2) sin(m phi_s)`.

The surface projection uses the outer-wall meridional normal/tangent. Because
the circumferential width is specified as a physical arc length, the raw
periodized-Gaussian angular integral is proportional to `sigma_s/r`; this
cancels the surface Jacobian `r` in `dS = r*ds*dphi`. The common factor
`sqrt(2*pi)*sigma_s` cancels in unit-impulse normalization, leaving the
meridional segment weight `ds` times the axial Gaussian. The normalized
Fourier moments above then supply the directional factors. No angular
quadrature is used in this projection. The resulting generalized impulse is
normalized to unit total scalar impulse and remains signed per eigenvector
convention.

The function rejects missing, malformed, non-finite, duplicated or unpaired
Fourier boundary data. It does not silently fall back to the axisymmetric
boundary arrays. The azimuthal displacement component is retained/validated
as part of the vector FEM record but does not contribute to a meridional
normal/tangent load because that load has no circumferential force component.

### Verified evidence

- Synthetic cosine/sine modal pair: rotation law and pair-norm invariance pass.
- `m=0` projection is invariant under strike azimuth.
- Increasing physical circumferential width attenuates the `m=1` contribution.
- Incomplete cosine/sine basis is rejected.
- The existing direct vector FEM evaluator successfully produces Fourier
  boundary data for an exactly circular geometry (`q=1`), and the new
  projection returns finite modal responses from those modes.
- On the tested `q=1` FEM mode window, changing only strike azimuth changes the
  modal generalized impulses and the reconstructed outer-wall response obeys
  the circular rotation covariance law.
- A three-mesh diagnostic (`nNodes=15,23,31`) gives relative outer-wall
  rotation-covariance errors between `2.6e-5` and `1.5e-4` for the retained
  16-mode window. The lowest `m=1` pair's generalized-force power changes by
  about `2.4%` between the coarsest and finest tested meshes. These are
  numerical observations, not acceptance thresholds.
- The same audit exposed and corrected a DOF-layout defect in the direct FEM
  evaluator and diagnostic modal state tracker. Assembly uses
  `[node][component][Fourier basis]`; the old evaluator read some fields as if
  the order were `[node][basis][component]`. This corrupted axisymmetric
  influence extraction, boundary Fourier exports and parts of axis constraints.
  The full-field tracker used a third, inconsistent indexing expression when
  matching different Fourier basis sizes. All now follow the assembly layout,
  with explicit regression checks.
- This direct vector FEM evaluation is currently an explicit API path. The
  default circular optimizer path still uses `InternalModalSolver` and its
  axisymmetric modes; the new strike projection is not wired into optimization
  or ordinary synthesis.

### Remaining gates

1. Complete contact-width and frequency convergence on a declared, sufficiently
   refined mesh sequence. The current mesh audit is diagnostic only: the lowest
   pair frequency still changes materially with refinement.
2. Extend covariance coverage to additional circular mode families, angular
   quadrature settings and eigenvector/subspace perturbations.
3. Revisit coordinate-axis regularity for vector Fourier modes. The current
   constraints now target the intended DOFs, but the physical regularity
   conditions for `m=0`, `m=1` and higher orders require a separate derivation.
4. Define impact impulse/time history and the resulting complex modal response;
   the current function returns only real signed generalized impulse factors.
5. Define and validate non-axisymmetric acoustic radiation coefficients.
6. Only then populate `AcousticDirectionalField` and implement stereo rendering.

This projection is a structural excitation step, not a radiation model and not
proof that the current synthesis path is directional. A separate normalization
contract test also verifies that fixed physical arc width cancels the local
circumference factor in the integrated patch area; equal meridional segment
lengths are therefore not spuriously weighted by radius.

## Explicit non-goals

This audit does not:

- change the production excitation;
- claim a physical striker width from existing diagnostic widths;
- reopen the isolated Chinese-bell research branch;
- infer acoustic radiation from structural Fourier coefficients;
- add stereo panning disguised as physical rotation.

## Status

**Localized strike projection:** implemented as an explicit direct-vector-FEM API; not wired into the default circular optimizer or ordinary synthesis.

**Non-axisymmetric structural excitation:** projection contract exists; circular rotation covariance is verified on actual FEM modes for the tested window. Contact-width/frequency convergence and complete axis regularity remain open.

**Directional acoustic field:** representation exists, population not implemented.

**Stereo:** blocked.

## 1.2.132: FEM Fourier DOF-layout correction and covariance audit

The direct vector FEM element/global assembly stores DOFs as `[node][component][Fourier basis]`. The evaluator had incorrectly interpreted several arrays as `[node][basis][component]`. This affected the generalized-influence vector, axisymmetric radial/axial traces, outer-surface Fourier export and coordinate-axis fixed-DOF indices. The full-volume modal state tracker also used an inconsistent basis-major index when comparing records with different Fourier basis counts.

The evaluator and tracker now use the assembly layout. The circular FEM integration test checks every exported outer-surface component/basis against `fem_displacement`, then reconstructs the modal outer-wall response for two strike azimuths and verifies rotation covariance. On `nNodes=15,23,31`, the relative covariance error was `1.44e-4`, `2.59e-5` and `8.10e-5`. The lowest `m=1` pair frequency moved from `5.13406 Hz` to `4.37053 Hz` (about `-14.9%`), while its summed generalized-force power changed about `-2.39%`. Therefore the rotational law is supported for this model/window, but the frequency mesh sequence is plainly not converged and is not a physical acceptance result.

A full-Fourier modal tracker regression also verifies that common components retain unit similarity when the candidate basis is enriched from `{0,1,-1}` to `{0,1,-1,2,-2}`.


## 1.2.133 axis regularity dependency

The circular strike projection depends on valid direct vector-FEM mode traces. The 1.2.133 fixed-DOF policy was incomplete, so its covariance numbers are provisional historical diagnostics only. Exact coupled-constraint elimination is implemented in the 1.2.134 follow-on; see the next section. This does not promote the projection to the production optimizer.


## 1.2.134 follow-on: exact axis regularity and retained mode window

The direct vector-FEM evaluator now uses exact `u=Pq` axis regularity elimination rather than individually pinning every non-axisymmetric trace. The `m=1` Cartesian smoothness constraints are `ur(1c)+uphi(1s)=0` and `ur(1s)-uphi(1c)=0`; this preserves the two constant transverse Cartesian traces. The previous 1.2.133 audit had the signs reversed in its displayed equations; that has been corrected.

The evaluator requests six free-body rigid modes plus one complete pair beyond the requested physical mode count, then applies the minimum-frequency filter. This prevents the output window from ending halfway through a degenerate circular pair. In the 15/23/31-node circular strike diagnostic, retained lowest-pair frequencies were 81.0672/79.7039/79.2442 Hz, rotation covariance error was 9.10e-5/5.62e-5/2.60e-5, and lowest-pair generalized-force power changed -32.6% from coarsest to finest. Frequency convergence improved to -2.25%, but the force-power sensitivity means contact/mesh convergence remains open.

The exact constraints do not populate `AcousticDirectionalField`, do not add impact time history, and do not make the direct vector-FEM path the default production synthesis path.
