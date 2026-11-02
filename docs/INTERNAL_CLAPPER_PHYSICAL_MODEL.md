# Internal clapper physical model and audio synthesis (1.2.140)

## Scope

`bell synth-clapper` is the first end-to-end internal-clapper path. It replaces the external translating-striker assumption for this command with a planar rigid physical pendulum whose spherical steel head contacts the true inner bell wall. The existing `synth-contact` command remains available as a separate external-striker diagnostic; it is not silently reinterpreted.

Pipeline:

1. Build the production bell mesh and solve direct cylindrical vector FEM with real Fourier components `m=0..3`.
2. Read the actual inner-wall modal displacement fields exported from the FEM eigenvectors.
3. Move a spherical head on a pendulum arc and find its closest point / signed gap against the SI inner meridian at every integration stage.
4. Evaluate inner-wall modal displacement at the contact segment and clapper swing azimuth.
5. Apply a unilateral Hertz force to the bell modes and the equal-and-opposite force to the pendulum rotation equation.
6. Stop at first separation, continue the bell modal state with per-mode T60 decay, and render the complex directional far-field radiation approximation from the outer wall to mono PCM16 WAV.

## Rigid pendulum

The model uses a planar swing plane at azimuth `phi_c`, with `z` positive from crown toward mouth. For pendulum angle `theta` measured from the downward vertical:

\[
r_h = r_p + L_h\sin\theta,\qquad z_h=z_p+L_h\cos\theta.
\]

The generalized coordinate is `theta`. The clapper mass properties are represented by total mass `m`, centre-of-mass distance `l_cm`, and pivot inertia `I_p`:

\[
I_p\ddot\theta=-m g l_{cm}\sin\theta + Q_c.
\]

The default inertia is the point-mass-equivalent estimate `I_p=m*l_cm^2`; it is only a starting estimate. For a real clapper, supply measured/calculated pivot inertia with `--clapper-inertia`. The head-centre offset `L_h` is independent of `l_cm`.

Initial angle and angular velocity are initial conditions, not inferred from bell geometry. The default example starts at the downward vertical with a positive angular velocity so the head swings toward the wall. A released-from-rest case can be specified with a nonzero initial angle, provided the direction and available gravitational energy lead to contact.

## Inner-wall geometry and modal coupling

The contact query uses `CylindricalBellModel::mesh_si.inner_profile`, not the outer wall and not a thickness approximation. The closest-point routine returns the cavity-facing normal `n_cav` and signed geometric gap `g = distance - R_tip`.

For each retained mode, the radial and axial Fourier fields are linearly interpolated on the contact segment and evaluated at `phi_c`. Their projection onto `n_cav` gives the local modal shape `phi_n`. With modal coordinates `q_k`, the wall displacement toward the cavity is:

\[
u_n=\sum_k \phi_{n,k}q_k,\qquad
\delta=\max(0,u_n-g).
\]

The unilateral elastic contact force is:

\[
F_c=k_H\delta^{3/2}.
\]

The modal generalized force is `Q_k = -phi_n,k F_c` because the head pushes the bell wall opposite the cavity-facing normal. The force on the head is `+F_c n_cav`; its pendulum generalized torque is `Q_c=F_c n_cav dot d(x_h)/d(theta)`. These signs make the work exchange between pendulum, bell modes and Hertz contact potential conservative in the undamped impact model.

## Hertz coefficient

The default coefficient is calculated from the steel head and the current bell material using the reduced elastic modulus:

\[
E^*=\left[\frac{1-\nu_h^2}{E_h}+\frac{1-\nu_b^2}{E_b}\right]^{-1},\qquad
k_H=\frac{4}{3}E^*\sqrt{R_h}.
\]

This is a **locally flat wall approximation**: it uses the sphere radius as the effective radius and does not yet include the two principal curvatures of the concave bell wall or a non-circular Hertz contact ellipse. The coefficient can be overridden for sensitivity studies with `--hertz-stiffness`, but it must not be treated as calibrated unless supported by measured contact data. The model is elastic; contact damping, plasticity, friction, wear and finite-element contact patches are not included.

## Acoustic output

At first separation, the computed modal displacement and velocity are continued with modal exponential decay using the requested T60. The outer wall remains the radiating surface. Each mode uses the existing complex directional far-field surface integral, preserving phase before summation. Output is peak-normalized mono PCM16 WAV. It is a repeatable directional listening preview, not calibrated SPL, BEM, room acoustics or stereo.

The short impact integration is conservative and undamped. Modal T60 decay is applied after separation, not during contact. The first version stops after the first complete contact/separation event; it does not simulate later rebounds or repeated strikes in the same file.

## Example

```text
bell synth-clapper bells/best_geometry.json --output audio/bell_clapper.wav \
  --clapper-mass 5.0 --com-length 0.46 --tip-length 0.52 --tip-radius 0.035 \
  --initial-angle-deg 0 --initial-angular-velocity 3.2 \
  --strike-azimuth-deg 0 --listener-polar-deg 30 --listener-azimuth-deg 0
```

All lengths are metres in the clapper command. `--pivot-z` defaults to 30 mm above the first inner-profile point; provide the actual suspension location for a measured bell. `--clapper-inertia` defaults to `m*l_cm^2` and should be replaced with a measured or CAD-derived value for calibrated work.

## Validation and limitations

- `clapper_dynamics_contract` checks that the pendulum reaches contact, separates, transfers energy to modal coordinates and conserves total undamped mechanical/contact energy within tolerance.
- `internal_clapper_audio_smoke` runs the geometry-to-WAV path on a coarse mesh.
- The bundled numerical parameter set (5 kg nominal clapper mass for this medium-sized example geometry) is a demonstrator, not a measured clapper specification; actual mass and inertia must be supplied for a particular instrument. The sound is not yet validated against a recording of a known bell/clapper pair.
- Mesh/modal spectral convergence, curvature-corrected Hertz contact, measured clapper mass properties, contact-patch convergence, repeated impacts and absolute acoustic calibration remain open.

## Numerical sweep on the bundled geometry

The coupled first-impact path was checked with the 5 kg demonstrator parameters:

- At 31 mesh nodes, reducing the integration step from 4 to 1 microsecond changed peak force by less than 0.002% (about 22.440 kN in all three runs). Contact duration remained 1.03–1.04 ms; conservative energy drift was about 3.15e-5.
- At 15, 31 and 71 mesh nodes with a 2 microsecond step, peak force was 22.775, 22.440 and 22.447 kN; contact duration was 1.01, 1.03 and 1.02 ms. First-contact time changed from 76.16 ms at 15 nodes to 75.66 ms at 71 nodes.

This is encouraging convergence of the **impact event metrics** for this parameter set. It does not establish modal spectral convergence, contact-patch accuracy, or acoustic waveform convergence; the structural eigenfrequencies and directional radiation still require their own mesh/family convergence audit.
