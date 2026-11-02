# Excitation models

## 1. Shared contract

The modal solver receives a generalized influence vector through `ModalInfluenceProvider`. The solver is not responsible for deciding whether that vector represents a strike, rotation proxy or another load case.

For an eigenmode `φ` and generalized force/influence vector `f`, the modal generalized force is

`Q = φᵀ f`.

If the mode is not mass-normalized, the participation normalization is `|Q| / sqrt(φᵀ M φ)`. A mass matrix must not be applied to the force vector before this projection.

## 2. Current optimizer influence

`src/excitation/meridional_influence.*` implements the existing meridional influence field used by the production optimizer. It is a proxy in the meridional `(r,z)` model; it is not a full three-dimensional impact or a claim of a complete rigid-body rotational excitation model.

## 3. Localized strike

`modalImpulse()` is the legacy meridional/ring-localized diagnostic. Its Gaussian depends on `z` and the surface measure integrates the full circumference, so it remains axisymmetric.

`modalImpulseCircular()` is a separate structural projection for a localized patch at `(z, phi)`. It requires complete direct-vector-FEM outer-surface Fourier data, uses a periodized Gaussian in physical circumferential arc length, and analytically projects that load onto cosine/sine basis labels. It returns signed real generalized impulse factors normalized to unit total scalar impulse. It does not model impact duration, striker compliance, complex modal phase, or acoustic radiation. The API is not yet wired into the default circular optimizer, whose current modal solve remains axisymmetric.

## 4. Chinese-bell strike

`research/chinese_bell/src/excitation/chinese_bell_strike.*` constructs the experimental contact distribution on the mapped outer surface. The force follows the local surface normal. For elliptic geometry the normal may have radial, azimuthal and axial components; the circular limit removes the azimuthal component.

The contact footprint uses segment-lumped surface quadrature. The spatial footprint is normalized as a force distribution; it does not model striker compliance, impact duration, contact stiffness or the energy transfer of a real mallet. Its width is an input in metres, and azimuthal width is in radians. Contact width should be grounded in measured or specified striker geometry, not optimized freely to maximize selectivity.

## 5. Acoustic interpretation

Modal participation is not radiated loudness. A complete evaluation must separately account for structural mode shape, generalized force, damping, acoustic radiation efficiency/directivity and observation position. The current scalar production path stores `Qd`, an axial first moment and a normalized pressure coefficient; it does not retain a spatial acoustic field. The current Chinese objective's monopole/dipole terms are diagnostic radiation proxies, not a calibrated non-axisymmetric acoustic solver.

The current scalar pressure coefficient combines the meridional participation `part_rz` with `Qd`. The meridional influence defaults to a `z_ref` of zero, while the acoustic first moment uses the bell mouth as its origin. Existing rotational-reference tests show that participation is origin-dependent. Therefore this coefficient must not be promoted to an origin-invariant physical pressure field without an explicit shared reference-frame contract.

For the planned stereo feature, the acoustic field must be represented with angular phase/directivity information before channel rendering. A circular-bell field expressed in azimuthal harmonics can then be rotated by phase factors `exp(-i m theta)` without another structural solve. See `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md`.

## 6. Production boundary

The Western production optimizer continues to use the meridional influence contract. The circular localized-strike projection is an explicit reusable API over direct vector-FEM output, not yet a default optimizer/synthesis path. Chinese strike and two-tone evaluation remain experimental and disconnected from the production optimizer until their independent validation gates pass.


## 1.2.132 FEM modal-record indexing correction

The direct vector-FEM global DOF layout is `[node][component][Fourier basis]`, not `[node][basis][component]`. The evaluator's generalized influence, axisymmetric traces, Fourier boundary export and axis-constraint indices now follow the assembly layout. The diagnostic `ModalStateTracker` full-volume comparison was corrected to use the same mapping when comparing common harmonics across enriched bases. Regression tests cover exported boundary data, enriched-basis tracking and actual circular strike-response rotation covariance.


## 7. 1.2.135 finite-duration modal response boundary

`modalImpulseCircular()` returns signed generalized impulse factors, not an
audio waveform. `modal_impulse_response.*` now provides the next independent
linear-dynamics layer. For each modal term it evaluates

`q_k(t) = J_k/(M_k*omega_d) * exp(-alpha_k*t) * sin(omega_d*t)`

with `alpha=ln(1000)/T60` and `omega_d=sqrt((2*pi*f_k)^2-alpha^2)`.
The implementation accepts only finite, positive parameters and underdamped
terms. It sums modes linearly, preserves the sign of the generalized impulse,
and validates the expected initial velocity `q'_k(0)=J_k/M_k`.

This is an ideal impulse response. It does not infer a real hammer contact
duration, contact stiffness, force envelope, nonlinear contact loss, or
frequency-dependent damping. Those remain separate measured/modelled inputs.
The current module is not wired into ordinary synthesis.

## 8. Contact and mesh convergence status

The 1.2.134 direct circular-vector-FEM diagnostic showed a relatively stable
lowest-pair frequency but a 32.6% coarsest-to-finest change in lowest-pair
generalized-force power. This is not an acceptable basis for calibration or
for claiming mesh-independent strike loudness. The diagnostic uses a fixed
physical strike definition; the remaining sensitivity may combine profile
quadrature, mesh-dependent modal traces, and eigenvector/subspace tracking.
It must not be hidden by renormalizing each mesh's modal response.

The older `localized_strike_mesh_convergence` test exercises the legacy
axisymmetric ring-localized API and is not evidence for circular point-strike
convergence. The actual circular integration test checks rotation covariance
on three meshes, but covariance alone is not convergence. The next numerical
gate is a controlled sweep that holds physical contact widths, strike location,
material and geometry fixed while refining the same circular FEM model, and
reports frequency, each degenerate pair's summed impulse power, and changes
under contact-width refinement. Until that evidence is available, physical
contact convergence remains **OPEN**.

## 9. Modal-family coverage

The projection contract now includes an explicit synthetic `m=2` cosine/sine
pair in addition to `m=1`. It verifies the `2*theta` phase law and pair-norm
invariance. This tests the Fourier convention at another order, but it does
not claim that an actual FEM `m=2` family has converged or that a full
multi-family physical response has been validated. Actual family tracking
and contact-power convergence remain required.

## 10. Finite striker contact (1.2.136)

`hammer_contact.*` adds a separate conservative point-contact integrator. The
striker is a free translational mass and the bell is represented by the supplied
linear modal coordinates. At the strike point,

`delta = x_striker - sum(phi_k(point) * q_k)`

and the unilateral Hertz force is `F = k_H * delta^(3/2)` for `delta > 0`,
otherwise zero. The contact potential `2/5 * k_H * delta^(5/2)` is included in
the energy audit. Equal-and-opposite interaction is applied: the striker
accelerates under `-F`, while mode `k` receives generalized force
`phi_k(point) * F`.

Inputs explicitly carry striker mass (kg), initial impact velocity (m/s), Hertz
stiffness (N/m^(3/2)), integration step (s), and modal frequency, modal mass,
and point mode-shape value. The solver uses fixed-step RK4. It reports striker
state, contact-force history, modal coordinates, total energy, maximum relative
energy drift and estimated contact duration.

This is a first conservative mechanical model, not a calibrated mallet model.
It omits contact damping, striker-tip layers, plasticity, finite contact-patch
geometry, friction, and any spatial variation of the contact point during the
impact. Hertz stiffness and modal point values must be physically calibrated;
modal mass and mode-shape normalization must be mutually consistent. The
integration step must be refined until force pulse, contact duration and energy
balance converge. The API is not wired into default synthesis, and the earlier
32.6% circular strike generalized-force mesh sensitivity remains OPEN.


## 11. Geometry-to-WAV contact synthesis (1.2.137)

The opt-in `bell synth-contact` command connects the Hertz contact state to captured axisymmetric FEM mode shapes and the existing outer-surface volume-displacement coefficients. It detects striker separation, continues the modal coordinates with a configured T60, and renders a monopole volume-acceleration proxy to mono PCM16 WAV. The solver conversion is explicit: the internal FEM uses mm and tonne-based mass normalization; contact and radiation values are converted to SI before integration/rendering.

This path is intentionally limited to the axisymmetric `m=0` projection of a point strike. It omits non-axisymmetric response and does not claim calibrated SPL or directional radiation. Contact metrics were stable across the tested 15–71 mesh-node sweep and 1–2 microsecond steps, but the modal spectral distribution remains mesh-sensitive. See `CONTACT_DRIVEN_AUDIO_PIPELINE.md`.


## 12. Fourier-localized contact and directional audio (1.2.138)

The opt-in `synth-contact` path now uses the direct cylindrical vector FEM rather than the axisymmetric internal solver. It solves the real Fourier basis `m=0..3`, exports the complete mass-normalized radial/azimuthal/axial surface fields, evaluates radial/axial displacement at the strike azimuth, and supplies that point-normal mode value to the conservative Hertz contact integrator. The modal mass and shape normalization are taken from the same SI direct FEM solve.

After separation, each mode is rendered through a complex far-field surface integral

`H_k(n) = integral_S (phi_k . n_surface) exp(-i k_k n . x) dS`.

The integral uses periodic azimuth quadrature and meridional segment integration. The analytic modal acceleration is multiplied by the complex transfer coefficient before taking its real part, preserving the frequency-dependent spatial phase at the selected listener direction. The output is still peak-normalized mono PCM16, not calibrated SPL.

A CLI CTest covers the Fourier FEM-to-WAV path. A manual rotational covariance check on the bundled circular geometry gave correlation 0.9999999985 when strike and listener azimuth were both rotated by 90 degrees; changing only the strike azimuth changed the waveform (correlation 0.803). This verifies a useful symmetry contract but does not establish physical acoustic accuracy. Harmonics are truncated at m=3, and the modal spectral mesh convergence / physical point-contact generalized-force sensitivity remain OPEN. Contact damping, finite tip geometry, validated radiation efficiency, stereo channel rendering, and room acoustics are not included.


## 13. Internal pendulum clapper coupling (1.2.140)

The new `synth-clapper` path is distinct from the external translational `synth-contact` path. `ClapperDynamicsConfig` carries clapper mass, pivot inertia, COM length, head-centre length, spherical-head radius, pivot coordinates, initial angular state and Hertz coefficient. `simulateInternalClapperContact` solves the planar physical-pendulum equation while querying the true SI inner profile and evaluating inner-wall modal fields at the moving contact point. The contact gap includes the modal wall displacement toward the cavity; the modal generalized force acts opposite the cavity normal, and the equal/opposite head force produces the pendulum torque. The impact is undamped and the total energy includes pendulum kinetic/gravitational energy, modal energy and Hertz potential.

The default Hertz coefficient uses `4/3 E* sqrt(R_tip)` with the reduced modulus of steel and bell material, assuming a locally flat wall. This deliberately remains an approximation: the concave inner wall's principal curvatures and elliptical contact patch are not yet accounted for. Modal T60 decay starts after first separation; only the first contact/release event is simulated. The generated file is peak-normalized mono directional far-field preview, not calibrated SPL. `clapper_dynamics_contract` and `internal_clapper_audio_smoke` cover the coupled model and audio path.
