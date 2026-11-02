# Validation and evidence

## Build contract

The supported release build is C++17 with no third-party C++ libraries. `build.cmd` configures and builds the project; `test.cmd`/`ctest` execute the registered tests. CalculiX is optional and is used only as an external differential reference.

## Validation layers

### Geometry and mesh

- profile generation and spline convergence;
- thickness and mesh convergence;
- castability and geometry-path contracts;
- input-safety and model-layer contracts.

### Axisymmetric modal solver

- eigenpair integrity;
- generalized residual and M-orthogonality checks;
- sparse solver scaling/robustness;
- modal objective convergence.

### Cylindrical vector FEM

- Q4 strain/assembly contracts;
- global modal solve and ellipticity convergence;
- cylindrical geometry path and reference-mesh contracts;
- cross-mesh/full-field modal comparison.

### Modal identity

- mass-weighted overlaps and principal-angle subspaces;
- conservative ambiguity handling;
- modal-family continuation;
- diagnostic state tracker with Kalman frequency prediction.

### Excitation and acoustics

- meridional influence;
- localized ring-strike direction/distribution and convergence;
- circular Fourier strike projection, cosine/sine rotation covariance, m=0 invariance, and q=1 direct-vector-FEM integration;
- direct-FEM modal-record DOF-layout checks, full-Fourier tracker basis-enrichment, and reconstructed outer-wall rotation covariance across three profile meshes;
- acoustic pressure scaling/phase and modal-pressure contracts;
- multipole and radiation diagnostics.

Passing these numerical contracts does not establish experimental calibration.

### External differential reference

CalculiX comparisons check implementation agreement on selected same-geometry cases. They do not establish historical or experimental accuracy.

## Current Chinese-bell status

The Chinese-bell implementation is diagnostic/experimental and remains disconnected from the production optimizer. The complete historical Chinese programme is preserved under `research/chinese_bell/`; start with `research/chinese_bell/docs/CHINESE_RESEARCH_HANDOFF_1.2.126.md`.

Durable conclusions:

- the corrected surface-integrated non-axisymmetric strike model is retained;
- modal identity must be based on spatial/subspace information, not frequency rank alone;
- earlier large frequency-split observations included modal-window/association artifacts;
- the current reduced FEM model has not established a physically identified `(m,n)a/(m,n)b` family with reciprocal A/B strike selectivity;
- contact width remains a mathematical scenario parameter;
- Chinese acoustic radiation/damping is not yet validated.

## Maintained Chinese numerical contracts

Only implementation-protecting tests remain in the maintained Chinese test surface:

- Chinese-bell research contracts are outside the normal validation suite. They are built only with `BUILD_CHINESE_RESEARCH=ON`; see `research/chinese_bell/CMakeLists.txt`.

The experiment-specific sweep matrices and branch-search executables were removed in 1.2.122 after their conclusions were consolidated into the Chinese audit document.

## Numerical evidence policy

A passing test establishes only the property it tests. Numerical PASS must not be interpreted as physical validation. Provisional Chinese questions remain blocked until independent reference or measurement evidence exists.

## Non-negotiable rules

- Do not tune gates to make a candidate pass.
- Do not weaken a gate solely because an optimization run fails it.
- Do not silently mix Western and Chinese objectives.
- Do not use an experimental evaluator as a production physical gate.
- Do not claim historical authenticity from geometric or acoustic heuristics alone.


## 1.2.132 focused release validation

The following 11 targeted CTest cases passed after the direct-FEM Fourier DOF-layout correction: `acoustic_directional_field_rotation_contract`, `localized_strike_direction_contract`, `circular_localized_strike_projection`, `direct_elliptical_fem_evaluation`, `localized_strike_mesh_convergence`, `cylindrical_global_assembly_contract`, `cylindrical_global_modes_contract`, `cylindrical_modal_family_tracking`, `cylindrical_bell_geometry_modal_tracking`, `cylindrical_modal_field_cross_mesh`, and `modal_state_tracker_contract`. The production `bell` executable also built successfully. The full all-target build/test suite was not completed in the available build window and is not claimed as passing.

The three-mesh circular strike diagnostic established rotation covariance for the tested model/window (relative errors `1.44e-4`, `2.59e-5`, `8.10e-5`). The lowest pair frequency changed by about `14.9%` from `nNodes=15` to `31`; this sequence is not frequency-converged. These are numerical implementation checks, not physical validation of striker dimensions or acoustic radiation.


## 1.2.134 exact axis regularity

The direct vector-FEM evaluator now uses exact `u=Pq` constraint elimination. `cylindrical_axis_regularity_contract`, `direct_elliptical_fem_evaluation`, `circular_localized_strike_projection`, and `modal_tracker_fem_continuation` pass. The 15/23/31-node circular diagnostic reports rotation covariance errors `9.10e-5`, `5.62e-5`, and `2.60e-5`; the lowest retained pair frequency changes by `-2.25%`, while its generalized-force power changes by `-32.6%`. This is improved numerical behavior, but contact/mesh sensitivity remains unresolved and the direct vector-FEM path remains diagnostic rather than the default production acoustic path.
