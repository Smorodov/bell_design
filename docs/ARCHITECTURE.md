# Architecture

## 1. Application flow

```text
CLI
 ├─ optimize ── BellOptimizer ── ModalSolver
 │                              ├─ internal sparse FEM
 │                              └─ optional CalculiX adapter
 ├─ synth ───── stored modal data ── BellSynthesis
 └─ synth-geometry ─ BellGeometry ── modal evaluation ── BellSynthesis
```

The direct cylindrical vector-FEM path is used by dedicated evaluation and diagnostic workflows and by the explicitly opt-in ellipticity evaluation bridge. It is not the production objective path for the Chinese-bell optimizer.

## 2. Domain and geometry

- `src/core/bell_types.h` defines shared geometry, mesh, material and modal data contracts.
- `src/core/bell_material.*` owns the single material-preset registry. Configuration selects a named preset; solver code receives only the resulting mechanical properties.
- `src/core/bell_generator.*` constructs and validates the production meridional profile/mesh.
- `src/models/bell_model_profile.*` supplies tradition-neutral and tradition-specific geometric priors.

Tradition profiles are input priors, not claims of historical authenticity.

## 3. Optimization and modal solving

- `src/core/bell_optimizer.*` owns the evolutionary search and the current Western tuning objective.
- `src/core/modal_solver.hpp` defines the solver interface and neutral excitation-provider contract.
- `src/fem/internal_modal_solver.*` implements the internal axisymmetric modal solver.
- `src/ccx/ccx_solver.*` implements the optional external CalculiX reference adapter.

The legacy reduced-order coordinate-morphing calculation remains in use by the existing optimizer and geometry-listening path. It is not equivalent to solving the morphed geometry with direct vector FEM. It must not be silently described as such.

## 4. Direct cylindrical vector FEM

`src/fem/` contains:

- `cylindrical_q4_element.*` — element kinematics and local stiffness/mass matrices;
- `cylindrical_global_assembly.*` — global sparse K/M assembly;
- `cylindrical_global_modes.*` — constrained generalized eigenproblem;
- `sparse_matrix.hpp`, `sparse_eigensolver.hpp` — internal sparse data structures and eigensolver;
- `cylindrical_bell_model.*`, `cylindrical_modal_evaluator.*` — geometry adapter and direct evaluation. The global vector-FEM DOF layout is `[node][component][Fourier basis]`, with component order `(ur, uphi, uz)`; exported full-volume and outer-boundary modal records preserve this exact order/basis mapping;
- `cylindrical_modal_tracking.*` — same-DOF-space mass-metric overlap, subspace comparison and exact nested-basis embedding.
- `cylindrical_modal_field.*` — independent full-field comparison across different Q4 meshes and nested Fourier bases; it is compiled into the test-only `bell_experimental` library, not `bell_core`. It evaluates each field in its own Q4 element, locates the corresponding physical meridional point in the other mesh, and integrates the density-weighted displacement product over the common domain and azimuth.

The FEM layer is tradition-neutral. Ellipticity enters geometry and assembly, not a post-hoc frequency correction.

## 5. Modal identity: three distinct responsibilities

These modules are intentionally not duplicates:

1. `core/modal_tracking.*` reorders legacy `ModalData` by available outer-profile shape features for the existing optimizer.
2. `fem/cylindrical_modal_tracking.*` compares vector-FEM modes/subspaces in a common DOF space using a supplied mass metric and supports exact nested Fourier embedding.
3. `core/modal_state_tracker.*`, built only into the test/diagnostic-only `bell_experimental` library, is a diagnostic two-track continuation aid. It predicts frequency with a constant-velocity Kalman model and associates candidate modes using shape/subspace evidence. It is not wired into the production optimizer and does not alter eigenpairs.

The state tracker itself still uses normalized profile-coordinate registration for boundary traces and row-sum lumped mass for same-mesh full-field features. Both the evaluator and state tracker use the global `[node][component][Fourier basis]` layout; cross-mesh full-field comparison remains the separate `cylindrical_modal_field` utility. The new `cylindrical_modal_field` utility provides a separate cross-mesh volume comparison using consistent density-weighted integration; it is not yet connected to the state tracker. Its current contract requires matching ellipticity and density and rejects insufficient common-domain coverage.

## 6. Excitation and acoustics

- `src/excitation/meridional_influence.*` provides the existing meridional generalized-influence proxy.
- `src/excitation/localized_strike.*` keeps the legacy ring-localized diagnostic and adds an explicit Fourier-based circular patch projection over direct vector-FEM modes.
- Chinese-bell strike/objective code is intentionally detached under `research/chinese_bell/` and is not compiled into the normal production targets.
- The research module consumes `bell_core` and may be enabled explicitly with `BUILD_CHINESE_RESEARCH=ON`.
- `src/synthesis/bell_synthesis.*` converts scalar modal data into audio using the current acoustic approximation.

The current production acoustic representation is scalar per mode: it retains modal radiation coefficients such as `Qd`, an axial first moment, and a normalized pressure coefficient. It is **not** a spatial acoustic field and therefore cannot by itself support physically meaningful stereo rotation. The new circular strike projection supplies structural generalized impulse factors only; it does not populate acoustic radiation coefficients.

The planned stereo architecture is documented in `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md`: first construct and retain a production directional acoustic field with explicit angular/phase information, then rotate that field in post-processing for left/right observations. An exactly axisymmetric (`m=0`) field is rotation-invariant, so meaningful stereo requires retained non-axisymmetric excitation/radiation information. Rotation must not trigger another FEM or acoustic solve. Near-field spatial acoustics is intentionally out of scope for the first implementation.

Mechanical generalized force, modal participation and radiated sound pressure are separate quantities. The current Chinese radiation proxy is not a calibrated full non-axisymmetric acoustic solution.

## 7. Validation-only code

`src/validation/nonaxisymmetric_reference.*` builds an independent 3-D reference mesh and is compiled only into the test-only `bell_validation` library. `src/ccx/` writes/runs/parses CalculiX jobs and remains in `bell_core` because the CLI exposes an optional CCX solver. Neither is part of the direct cylindrical production FEM formulation.

## 8. Dependencies and build boundaries

The application and internal solver require only C++17 and the standard library. CalculiX is optional. `BUILD_SMOKE_TESTS=OFF` excludes all validation executables and CTest registration; it does not remove runtime application code.

## 8. Direct vector-FEM DOF layout contract

The direct cylindrical Q4 element orders each node's DOFs as `[ur basis 0..B-1][uphi basis 0..B-1][uz basis 0..B-1]`. Thus the global index is `node * (3*B) + component * B + basis`. `ModalMode::fem_displacement` and `fem_mass_lumped_weights` preserve this layout; outer Fourier traces use `[profile point][basis]` per component. The evaluator, axis-constraint construction and full-field modal tracker must use this same mapping.


## Cylindrical vector-Fourier axis regularity (1.2.134)

Smooth Cartesian vector fields require `m=0` transverse traces to vanish; `m=1` transverse traces to satisfy `ur(1c)+uphi(1s)=0` and `ur(1s)-uphi(1c)=0`; `m=1` axial traces to vanish; and all `m>=2` component traces to vanish. The direct vector-FEM evaluator now imposes these conditions by exact `u=Pq` elimination and forms `P^T K P` / `P^T M P`, preserving symmetry without penalty stiffness. See `docs/CYLINDRICAL_VECTOR_AXIS_REGULARITY_AUDIT.md` for the derivation and regression coverage.


## Contact-driven audio output (1.2.138)

The opt-in `synth-contact` path uses the direct cylindrical vector-FEM modal solve with a bounded real Fourier basis m=0..3. It evaluates each retained mode at the strike azimuth, couples the resulting point-normal values to the finite Hertz hammer model, and continues the separation state with modal T60 decay. A separate acoustic step integrates the complex far-field surface transfer coefficient for one listener direction; the complex phase is retained when converting modal analytic acceleration to a real mono signal. The result is peak-normalized PCM16 WAV, not calibrated SPL. The legacy/default synthesis path remains unchanged. See `CONTACT_DRIVEN_AUDIO_PIPELINE.md` and `EXCITATION_AUDIT.md`.

Structural Fourier displacement remains distinct from the acoustic transfer coefficient: the latter is computed by a spatially phased surface integral and is not stored as if it were a modal boundary field.


## Internal clapper mechanics and audio path (1.2.140)

The direct cylindrical modal evaluator exports mass-normalized Fourier displacement components on both the outer and inner meridional boundaries. Inner-wall arrays are copied from the actual inner-wall nodes in the full FEM eigenvector and are ordered like `inner_profile`; they are not inferred from the outer boundary. `Excitation::find_spherical_tip_contact_geometry` locates the spherical head against the cavity-facing inner meridian. In 1.2.140, `Excitation::simulateInternalClapperContact` adds a planar rigid-pendulum state and equal-and-opposite conservative Hertz coupling to the inner-wall modal coordinates. The new `synth-clapper` command sends the post-separation modal state through the outer-wall directional radiation approximation. The older `synth-contact` command remains a separate external translating-striker diagnostic. See `INTERNAL_CLAPPER_PHYSICAL_MODEL.md` for assumptions and limitations.
