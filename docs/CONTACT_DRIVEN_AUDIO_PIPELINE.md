# Contact-driven audio pipeline

## Scope

`bell synth-contact` is the first complete geometry-to-WAV path driven by a finite hammer impact. It is a physically motivated listening model, not a calibrated prediction of sound pressure level or a full three-dimensional point-contact solution.

Pipeline:

1. Load and validate a circular bell geometry.
2. Build the meridional Q4 mesh and solve the direct cylindrical vector-FEM problem with the real Fourier basis `m=0..3` and exact axis regularity constraints.
3. Evaluate every retained modal surface field at the selected strike height and azimuth, then project its radial/axial displacement onto the local surface normal.
4. Use the direct solver's SI mass-normalized modal coordinates consistently in the Hertz contact system.
5. Integrate the unilateral Hertz contact between a free striker and the retained Fourier modal system.
6. Detect separation and continue each modal coordinate as an underdamped free vibration with the requested T60.
7. Integrate the complex far-field surface-radiation coefficient for the selected listener polar and azimuth angles, including spatial phase delay; combine it with each modal analytic acceleration.
8. Peak-normalize and write mono PCM16 WAV plus a text metadata sidecar.

## Command

```text
bell synth-contact [geometry.json] [--output file.wav]
  [--strike-height 0..1] [--t60 seconds] [--duration seconds]
  [--sample-rate hz] [--mesh-nodes n]
  [--strike-azimuth-deg deg]
  [--listener-polar-deg deg] [--listener-azimuth-deg deg]
  [--striker-mass kg] [--striker-velocity m/s]
  [--hertz-stiffness N/m^1.5] [--contact-step seconds]
  [--contact-duration seconds]
```

Windows convenience script: `generate_contact_audio.cmd`.

## Physical and numerical conventions

- Contact law: `F = k * max(delta, 0)^(3/2)`.
- During contact, the striker and retained undamped linear modes exchange energy. The contact integrator records modal displacement and velocity so the free response starts from the actual separation state.
- Modal decay is applied after separation using `alpha = ln(1000)/T60`.
- The pressure proxy is proportional to the sum of each mode's surface volume-displacement coefficient times modal acceleration. The final digital peak is set to 0.95; this is not an SPL calibration.
- The default geometry is circular. Non-circular geometry is rejected rather than silently passed through an axisymmetric solver.
- Structural excitation uses the retained real Fourier basis `m=0..3`. Each mode's point displacement is evaluated at the specified strike azimuth; the hammer therefore couples to non-axisymmetric mode shapes rather than only the ring-average component.
- Radiation is a complex far-field surface integral for one observer direction. The surface integral uses a linear acoustic phase factor `exp(-i k n·x)` and a 96-point periodic azimuth quadrature. It is a first-order far-field approximation, not a boundary-element solution or calibrated SPL.
- The listener polar angle is measured from the bell axis (0° is above/below the axis, 90° is in the radial plane); listener azimuth is measured around the bell axis.
- The command remains restricted to circular geometry. Elliptical geometry requires a separate treatment of the azimuthally varying meridional shape and surface normal.

## Important physical limitation: not an internal clapper

The current command still represents a freely translating point striker and evaluates contact on the outer wall. It is not the suspended internal tongue/clapper used in a traditional bell. Its `--striker-mass`, `--striker-velocity`, and `--hertz-stiffness` are diagnostic parameters, not a calibrated clapper description. Do not interpret this path as physically correct for an internal spherical steel head. The required replacement is specified in `docs/CLAPPER_GEOMETRY_AND_KINEMATICS_AUDIT.md`. In particular, contact must be detected against the actual normal-offset `inner_profile`, use inner-wall modal displacement, and couple the equal-and-opposite contact force to rigid-pendulum motion. The outer wall remains the acoustic radiation surface.

## Validation status and limitations

The conservative contact integrator has a dedicated contract test for energy balance, striker impulse balance, invalid inputs, and time-step refinement. A CLI integration smoke test now exercises the full Fourier FEM-to-WAV path. On the bundled circular geometry, simultaneous rotation of strike and listener azimuth by 90° produced waveforms with correlation 0.9999999985; rotating only the strike while holding the listener fixed changed the waveform (correlation 0.803). A default 71-node run completed in about 10 seconds in the development environment. These are implementation/covariance checks, not physical validation. Modal spectral mesh convergence, contact-width sensitivity, and measured parameter calibration remain open.

Not included: contact damping, plasticity, finite striker tip geometry, material calibration from measurements, harmonics above m=3, validated absolute radiation efficiency, stereo channel rendering, room acoustics, and calibrated absolute sound pressure. The generated file is a repeatable model-based listening preview. Do not interpret its normalized amplitude as the loudness of a real bell.
