# Architecture

## Runtime paths

```text
CLI
 ├─ optimize ─ BellOptimizer ─ ModalSolver
 │                         ├─ internal axisymmetric solver
 │                         └─ optional CalculiX adapter
 ├─ synth ─ stored modal data ─ BellSynthesis
 ├─ synth-geometry ─ profile/mesh ─ modal evaluation ─ BellSynthesis
 ├─ synth-contact ─ direct cylindrical modes ─ striker/contact ─ modal release
 └─ synth-clapper ─ direct cylindrical modes ─ pendulum/contact ─ modal release
                         └─ shared far-field directional-field builder
                            ├─ mono receiver or left ear
                            └─ right ear
 All synthesis commands ─ shared PCM16 WAV writer
```

The production optimizer and the direct cylindrical vector-FEM listening path are distinct. The latter must not be described as a 3-D FEM solution or as the optimizer's default solver.

## Module ownership

- `src/core/` owns shared domain types, material presets, profile generation, optimizer and legacy modal tracking.
- `src/models/` owns geometry priors and bell-profile representations.
- `src/fem/` owns Q4 element formulation, cylindrical global assembly, sparse eigensolver, direct vector-FEM modal evaluation and same-space modal comparison.
- `src/math/` owns spline, coordinate-morphing and strain utilities. Coordinate morphing remains an active reduced-order path; it is not equivalent to solving the transformed geometry with direct FEM.
- `src/excitation/` owns localized excitation, external striker contact, inner-clapper geometry, pendulum dynamics, modal coupling and modal decay.
- `src/acoustics/` owns the direct far-field azimuthal reference integral, the production Fourier/Bessel directional-field builder, listener/ear geometry and finite-distance free-space receiver integral. The direct integral is retained as an independent numerical reference, not used by production synthesis.
- `src/audio/` owns the shared little-endian PCM16 RIFF/WAVE serializer for mono and stereo outputs.
- `src/synthesis/` owns modal synthesis; `src/commands/` owns CLI orchestration.
- `src/ccx/` owns the optional CalculiX adapter; `src/validation/` owns the independent reference-mesh path.

## Modal identity is not one algorithm

Three responsibilities remain separate:

1. `core/modal_tracking.*` tracks legacy `ModalData` using profile features for the existing optimizer.
2. `fem/cylindrical_modal_tracking.*` compares direct-FEM modes/subspaces in a common DOF space using a mass metric and exact nested Fourier-basis embedding.
3. `core/modal_state_tracker.*` is a diagnostic two-track predictor in the `bell_experimental` test-only library. It is not connected to the production optimizer. Full Q4 cross-mesh field comparison (`fem/cylindrical_modal_field.*`) is a separate diagnostic utility because the state tracker does not own Q4 topology.

The direct vector-FEM global DOF layout is `[node][component][Fourier basis]`, with component order `(ur, uphi, uz)`. Any exported volume or boundary modal records must preserve this mapping.

## Acoustic and stereo ownership

The `synth-clapper --stereo` path computes one structural solution, one contact event and one modal release history. Both channels evaluate that shared history at separate receiver positions. No FEM/contact/time-integration solve is performed per channel.

- `far-field` builds one complex azimuthal field per mode through `buildDirectionalRadiationField` and evaluates it at receiver directions; this is the default low-cost path. The mono `synth-contact` path and stereo/clapper far-field path now share this implementation.
- `finite-distance` evaluates an approximate free-space Green-kernel surface integral at each ear; it is an opt-in comparison path, not a full exterior acoustic solution.
- The receiver frame includes head yaw, ear spacing, ear positions and propagation offsets. Head diffraction, pinna/HRTF, torso scattering, room response and calibrated SPL are not implemented.

A rotation of the listener frame must not change the bell's FEM solution, clapper contact history or modal decay. Only receiver evaluation changes.

## Production and research boundaries

The Chinese-bell implementation lives under `research/chinese_bell/` and is built only with `BUILD_CHINESE_RESEARCH=ON`. It must not leak into the Western optimizer objective. Its current physical modal-family identification is unresolved; see its status document.

The clapper acoustic objective is a tested primitive but is not wired into production candidate evaluation. Any such integration requires a separate explicit design decision and physical calibration gates.

## Build dependencies

The normal application requires C++17 and the standard library. CalculiX is optional. `BUILD_SMOKE_TESTS=OFF` excludes test/diagnostic libraries and test executables. The Chinese research target is separately controlled and defaults OFF.
