# Clapper geometry and kinematics audit

## Status

Design audit for the next contact-model revision. This document does not change the production solver or claim that the current `synth-contact` output represents a real bell clapper. The current 1.2.138 contact path remains a free translational striker with a point Hertz contact.

## 1. Physical system to represent

The intended excitation is a suspended bell clapper (tongue) with a rounded, approximately spherical soft-steel head moving inside the bell. The head is carried by a pendular body and contacts the inner wall. It is not an externally approaching hammer.

The minimum useful model has three coupled parts:

1. a rigid-body clapper with a pivot, mass distribution and head-centre offset;
2. unilateral contact between the spherical head and the bell's inner surface;
3. the bell's elastic modal response, with contact force applied to the inner boundary.

The model should initially remain planar in one swing plane, with the swing-plane azimuth as an input. A freely moving 3-D clapper, flexible suspension, friction, wear and plasticity are later extensions, not prerequisites for a correct first model.

## 2. Existing geometry conventions and findings

`BellGeometry` uses millimetres. `z=0` is the crown/axis end and `z=H` is the mouth; the meridional profile is represented by `(r,z)`. `BellGenerator` constructs both `outer_profile` and `inner_profile`.

The inner profile is a normal offset of the outer meridional curve. Its points generally have both a changed radius and a changed axial coordinate. Therefore:

- collision detection must use the actual `inner_profile` polyline (or a documented smooth interpolation of it);
- it must not treat `getInnerRadius(z)` as the exact inner boundary, because that helper returns a radial estimate at the same `z`, whereas the stored inner profile is a normal offset;
- the contact normal must be computed from the local tangent of `inner_profile`, with orientation into the bell cavity;
- a normal computed from `outer_profile` is not a substitute for the inner-wall contact normal, especially on a coarse mesh.

The production mesh stores outer boundary nodes first and inner boundary nodes second. The direct cylindrical FEM evaluator already retains the complete mass-normalized eigenvector in `ModalMode::fem_displacement`, with DOF layout `[node][component][Fourier basis]` and component order `[ur, uphi, uz]`. It currently exports Fourier boundary arrays only for the outer surface. Thus the inner-wall modal data can be exported from the existing eigenvector and mesh ordering without another FEM solve; this is a focused evaluator extension, not a new structural model.

The current `synth-contact` command selects a point from `outer_profile`, evaluates `outer_surface_fourier_*`, and projects onto an outer-surface normal. That is explicitly the wrong boundary for an internal clapper impact. The existing directional radiation integral, by contrast, should continue to use the outer radiating surface.

## 3. Coordinate and contact convention

Use a right-handed cylindrical frame with `z` positive from crown toward mouth. Let the clapper swing in a fixed vertical plane of azimuth `phi_c`.

At a candidate inner-wall contact point `s`, define:

- `x_s(s, phi_c)`: point on the revolved inner profile;
- `n_cav(s, phi_c)`: unit normal from the bell material into the cavity;
- `x_h(theta)`: spherical head centre for clapper angle `theta`;
- `R_h`: spherical head radius;
- `g(theta,s) = |x_h(theta)-x_s(s)| - R_h`: undeformed geometric clearance at the closest point.

At smooth contact, the closest-point condition is normal to the inner profile and `x_h - x_s` is parallel to `n_cav`. Contact is active only when the minimum clearance reaches zero. The head must remain on the cavity side of the wall; a nearest-point search alone is not sufficient if it selects the wrong side of a folded or invalid profile.

For the modal contact convention, let `phi_n = phi_mode · n_cav` be the inner-surface modal displacement projected toward the cavity. A positive wall displacement in this direction moves the wall toward the head and increases compression. With positive compression

`delta = max(0, -g + u_n)`,

where `u_n` is the modal wall displacement toward the cavity, the force exerted by the head on the bell points opposite `n_cav`. The generalized modal force therefore carries the corresponding negative sign. The equal-and-opposite force on the clapper must be used in the rigid-body equation. This sign convention must be tested by total-energy conservation before production use.

The sign must not be hidden by taking an absolute value of the modal projection. Eigenvector sign is arbitrary, but the paired contact kinematics and generalized force must remain internally consistent.

## 4. Clapper rigid-body model

A useful first model is a rigid pendulum, not a point mass unless the clapper mass distribution is genuinely negligible.

Required physical parameters:

- total clapper mass `m_c`;
- pivot-to-centre-of-mass distance `l_cm`;
- pivot moment of inertia `I_p` about the swing axis;
- pivot-to-spherical-head-centre distance `L_h`;
- spherical head radius `R_h`;
- pivot position `(r_p,z_p)`;
- swing-plane azimuth `phi_c`;
- release angle `theta_0` or initial angular velocity;
- optional pivot damping, initially disabled unless measured.

For a simple pendulum angle measured from the downward vertical, the gravity torque is `-m_c g l_cm sin(theta)`. During contact, the contact force contributes a torque about the pivot. The equation is

`I_p theta_ddot = -m_c g l_cm sin(theta) + tau_contact`.

The head-centre position is obtained by rotating the pivot-to-head vector with the same rigid-body angle. Do not assume `L_h = l_cm`; they are different geometric quantities in a real clapper.

If released from rest at `theta_0`, the pre-contact angular speed in the lossless free swing follows energy conservation:

`0.5 I_p theta_dot^2 = m_c g l_cm (cos(theta) - cos(theta_0))`.

This relation is a useful independent check, not a replacement for time integration through contact. If the clapper is initially given a measured angular speed, use that directly instead.

The clapper is not automatically a point mass, nor is its impact speed an independent free parameter once pivot geometry and release conditions are specified. A compatibility mode may later accept measured impact speed, but the physical pendulum path should derive it from the rigid-body state.

## 5. Contact law and material model

The existing law `F = k_H delta^(3/2)` is a reasonable first elastic local-contact approximation, but its coefficient must be derived from the head radius, signed local curvature of the inner bell surface, and both materials. It must not remain an unexplained global tuning constant in the physical clapper mode.

For a sphere against a locally flat elastic half-space, the classical coefficient is

`k_H = (4/3) E_star sqrt(R_h)`,

with

`1/E_star = (1-nu_h^2)/E_h + (1-nu_b^2)/E_b`.

The bell's inner surface is curved and concave. The actual local contact is generally elliptical; its relative curvature must be formed from the signed principal curvatures of the sphere and inner bell surface in both tangent directions. The flat-surface expression is only a limiting check. Curvature signs must be handled geometrically, not guessed from unsigned radii.

Hertz contact is only appropriate while deformation is elastic, the contact patch is small relative to local curvature radii and wall thickness, and half-space assumptions are adequate. If the contact patch becomes comparable to wall thickness or local curvature scales, shell compliance and a finite contact patch need a better model. A soft-steel head is not a rigid body at the contact interface merely because its Young's modulus is high; both material elastic properties enter the reduced modulus. Plasticity and local damage should not be added without evidence that the elastic regime is exceeded.

The first implementation should expose `R_h`, `E_h`, `nu_h`, and the bell material properties explicitly, calculate the contact coefficient (or documented elliptic-contact equivalent), and report the assumptions used. No numerical defaults should be presented as measured traditional-clapper values.

## 6. Coupling to the bell modes

For each mode and candidate inner contact point:

1. interpolate the inner-boundary Fourier displacement from the full FEM eigenvector;
2. evaluate the real Fourier field at the clapper swing azimuth;
3. project `(ur, uz)` onto `n_cav` (the azimuthal component is zero for a purely meridional normal, but it should remain available for a future non-axisymmetric contact surface);
4. use that signed modal projection in the contact gap and generalized force;
5. apply the equal-and-opposite contact force to the clapper's rotational equation.

The contact point is not fixed at a user-selected normalized outer-profile height. It is found from the head-centre trajectory and the inner surface. At first contact, the nearest point and its normal define the local contact frame. During a small contact interval, the contact point may move along the surface; the first model may use a consistent closest-point update each integration step, provided it is stable and its energy balance is demonstrated.

The FEM mode shape used for contact must be sampled on the inner surface. The outer-surface mode data remains appropriate for the current outer radiation integral. These two surfaces must not be conflated.

## 7. Required numerical contracts

Before generating production audio, add independent tests for:

- **Geometry:** a sphere at known clearance from a simple circular/analytic inner profile; correct closest point, cavity-side test and contact normal.
- **Rigid-body kinematics:** head-centre trajectory agrees with the rigid transform; zero-angle equilibrium and free-pendulum energy are correct.
- **Impact timing:** first contact is detected at the geometric event, not at an arbitrary preselected time step.
- **Action/reaction:** clapper angular momentum change and bell generalized impulse have consistent signs and magnitudes.
- **Energy:** in the undamped elastic system, clapper kinetic/potential energy + bell modal energy + Hertz contact potential remain conserved to a stated tolerance.
- **Time-step convergence:** contact time, peak force, impulse and modal energy converge independently of the audio sample rate.
- **Mesh convergence:** inner-wall contact projection and modal frequencies converge under meridional refinement; compare degenerate modal subspaces, not only mode indices.
- **Rotation covariance:** rotating the swing-plane azimuth and all observer azimuths together on a circular bell leaves the appropriately rotated result invariant.

The current 1.2.138 energy test for a translating free striker is not sufficient evidence for a pendulum contact implementation.

## 8. Recommended implementation order

1. Export inner-surface Fourier modal displacement arrays from the existing full FEM eigenvector; add a boundary-index contract test.
2. Add a pure geometry utility for closest point, signed gap and cavity-oriented normal on `inner_profile`; test it without FEM.
3. Add a rigid-clapper kinematics type and free-pendulum tests; do not yet couple it to the bell.
4. Add coupled clapper/modal contact with equal-and-opposite forces and a total-energy audit.
5. Derive the Hertz coefficient from head and bell parameters; add curvature/validity diagnostics.
6. Only after the above pass, route `synth-contact` through the clapper model and label the old translating-striker path as a separate diagnostic if it is still useful.

## 9. Inputs still required for a calibrated example

No single set of traditional clapper dimensions should be invented. A calibrated case needs at least:

- clapper total mass and either its centre of mass plus pivot inertia, or its full mass geometry;
- pivot location and pivot-to-head-centre distance;
- spherical head radius;
- release angle or measured angular/impact velocity;
- steel-head elastic constants and bell alloy elastic constants;
- confirmation of whether the bell is stationary and the clapper swings in a fixed plane, or whether the bell itself swings.

If these are unknown, the implementation can still support parameter sweeps and sensitivity analysis, but output must be identified as a parametric model rather than a prediction of a specific historical bell.

## 10. Explicit non-goals for the first clapper revision

Do not add 3-D FEM, flexible clapper-body modes, contact friction, wear, plasticity, nonlinear bell material, room acoustics or stereo as part of this step. Do not change the detached Chinese-bell research tree. Keep the generic structural solver independent of clapper details; the clapper belongs in the excitation/contact layer.

## 11. Implemented foundation in 1.2.139

Two preparatory contracts from the audit are now implemented without changing the production `synth-contact` path:

1. `ModalMode` carries `inner_surface_radial`, `inner_surface_axial`, and full Fourier arrays `inner_surface_fourier_{radial,azimuthal,axial}`. The direct cylindrical modal evaluator copies these values from the actual inner-wall nodes in `fem_displacement`, preserving the solver's mass normalization and harmonic ordering. The arrays follow `inner_profile` order; the scalar radial/axial arrays are the zero-harmonic view.
2. `Excitation::find_spherical_tip_contact_geometry` finds the closest point on the piecewise-linear inner meridian to a spherical tip centre, returns the segment and interpolation parameter, the cavity-facing normal, centre distance, and signed gap (`distance - tip_radius`). The query rejects a centre on the exterior side of the locally oriented inner wall and rejects invalid/degenerate inputs.

The geometric function is intentionally a meridional, axisymmetric nearest-profile calculation. It does not yet calculate a 3-D contact patch, curvature-corrected Hertz coefficient, elastic wall displacement, contact force, friction, or pendulum motion. At profile corners and endpoints, the nearest-polyline result is a geometric approximation whose normal should be refined with a smooth local profile representation before high-accuracy contact dynamics.

New contracts:

- `inner_surface_modal_export_contract`: checks every exported Fourier component against the corresponding full FEM eigenvector DOF at the inner-wall nodes.
- `internal_clapper_geometry_contract`: checks contact projection, cavity-normal orientation, positive/negative gap, exterior-side rejection, and invalid inputs.

The old 1.2.138 `synth-contact` remains unchanged and is still an external translating-striker diagnostic. It must not be described as an internal clapper simulation.


## 12. Coupled implementation in 1.2.140

The first production excitation implementation is now `Excitation::simulateInternalClapperContact` in `src/excitation/clapper_dynamics.*`. It uses the explicit physical-pendulum state `(theta, theta_dot)`, actual inner-wall modal Fourier fields at the current closest-point segment and swing azimuth, unilateral Hertz force, and equal/opposite generalized torque/modal force. The impact integrator reports first contact, first separation, peak force, modal release state and conservative energy drift. `bell synth-clapper` then continues those modal coordinates with T60 decay and renders through the outer-wall directional far-field approximation.

The default Hertz coefficient is computed from steel and bell reduced elastic modulus with `R_eff = R_tip`, i.e. a locally flat wall approximation. It is not curvature-corrected and must not be represented as a validated contact law for a concave shell. Default mass, inertia estimate, tip dimensions, pivot offset and angular velocity are illustrative. The first implementation stops at first separation and does not simulate subsequent rebounds.

New checks: `clapper_dynamics_contract` verifies contact/release, modal transfer and total conservative energy drift; `internal_clapper_audio_smoke` exercises the end-to-end geometry-to-WAV path. See `docs/INTERNAL_CLAPPER_PHYSICAL_MODEL.md` and `RELEASE_NOTES_1.2.140.md`.
