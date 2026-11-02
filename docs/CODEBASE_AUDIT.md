# Codebase audit — baseline 1.2.88, extended in 1.2.89

## Consolidation decisions

- Replaced the chronological README with a description of the current application, module map, production boundary and build workflow.
- Rewrote architecture, mathematical foundations, validation, development plan, excitation audit and user guide against the current source tree.
- Consolidated version-specific implementation notes into the current architecture, validation and Chinese-bell status documents. Removed superseded version-stamped documents rather than presenting them as parallel current specifications.
- Merged the two nearly identical mesh-convergence command wrappers into `test_convergence.cmd [all|mesh|thickness]`.
- Removed generated CSVs from the source/configuration directory and removed the old `audio/verified_1.2.39` configuration snapshot. The circular listening WAV and its generation metadata remain.
- Corrected the CLI version string to match CMake/project version `1.2.88`.

## Production/experimental boundary

`bell_core` contains the application and maintained physical/solver implementation. The diagnostic-only `ModalStateTracker` remains in `bell_experimental`. The Chinese-bell experimental implementation has now been detached into `research/chinese_bell/` and is not part of `bell_core` or the normal validation build. The independent 3-D reference mesh builder was moved into `bell_validation`, also test-only. The production executable does not link either library. A symbol check of `libbell_core.a` confirms that neither diagnostic component nor the reference mesh builder is present in the production library.

The following similarly named modules were deliberately retained because their responsibilities and call sites differ:

- `core/modal_tracking`: legacy `ModalData` ordering used by the current optimizer;
- `fem/cylindrical_modal_tracking`: mass-metric mode/subspace comparison for the vector FEM path;
- `core/modal_state_tracker`: diagnostic Kalman-assisted continuation;
- `math/coordinate_morphing` and `math/modal_coordinate_morphing`: still called by the existing optimizer fallback and geometry-listening command;
- `excitation/meridional_influence`: still supplies the production optimizer's current generalized influence;
- `ccx` and `validation`: optional independent reference path, still used by dedicated tests and the solver interface.

No active production function was deleted merely because its name resembled another API. Removing coordinate morphing or the CalculiX adapter now would break live call paths or remove an active independent validation route. Those removals remain gated by migration/validation in `DEVELOPMENT_PLAN.md`.

## Build-system corrections

- All validation executables and all `add_test` registrations now live inside one `BUILD_SMOKE_TESTS` block.
- `BUILD_SMOKE_TESTS=OFF` builds only `bell_core` and `bell`; it does not build `bell_experimental` or test executables.
- Removed the broad `-Wno-unused-parameter` suppression. The clean GCC build produced no warnings with `-Wall -Wextra -Wpedantic`.
- The CLI and CMake project version both report `1.2.88` in the 1.2.88 baseline.

## Validation evidence

- Full application/core build: PASS.
- Application-only build with `BUILD_SMOKE_TESTS=OFF`: PASS.
- CLI `--version`: `BellOptimizer 1.2.88`.
- The 1.2.88 baseline registered 68 tests. Its evidence is retained as historical release evidence.
- GCC warning scan after clean rebuild/resume: 0 warnings.

## Remaining intentional limitations

The Chinese two-tone objective is not connected to production optimization. Modal continuation across different meshes now has an independent full-volume comparison utility in 1.2.89, but the utility is not yet connected to the state tracker. The state tracker still uses its previous same-layout/lumped-mass path. The reduced-order coordinate-morphing path remains active until its callers are migrated. CalculiX remains an optional external differential-validation adapter.


## 1.2.89 cross-mesh modal field comparison

- Added `src/fem/cylindrical_modal_field.*` to the test-only `bell_experimental` library as an independent diagnostic numerical utility; it is not wired into `ModalStateTracker` or the optimizer.
- Added `cylindrical_modal_field_cross_mesh`, covering different Q4 meshes, nested Fourier bases, two-dimensional subspaces, common-domain rejection and actual FEM eigenmodes.
- Targeted cross-mesh test result: analytic overlap 1.0, two-dimensional subspace similarity 1.0, actual FEM best overlap 1.0; a clipped domain is rejected at common-volume fraction about 0.8611.
- The current comparison requires equal ellipticity and density. It uses 3x3 meridional Gauss quadrature and azimuthal trapezoidal quadrature; its continuum mass integral is not bitwise identical to the assembled FEM mass matrix (2x2 meridional element quadrature).
- CMake/CLI version is 1.2.89. Both full build and `BUILD_SMOKE_TESTS=OFF` application-only build pass.
- CTest registers 69 tests. 68 have been individually or batch-confirmed passing. `cylindrical_unified_ellipticity_audit` (#59) exceeded the available execution window and remains unconfirmed in this run; therefore the full suite is not declared green.


## 1.2.92 real-FEM modal continuation matrix

Added `cylindrical_modal_field_fem_matrix` to exercise the 11th/12th nearly-degenerate family from the Chinese-bell readiness study. It compares the full-field two-dimensional subspace across (a) mesh refinement only, (b) Fourier enrichment only, and (c) both changes simultaneously. The measured principal overlaps are 0.981226, 0.999712 and 0.980768; common-domain coverage is 1.0 in all three cases. The legacy boundary matcher reports no shape assignment in the mesh-only case. The new utility remains isolated in `bell_experimental`; no production tracker or optimizer behavior changes in 1.2.92.


## 1.2.91 memory and lifetime hardening

The 1.2.92 failure report localized two Windows `0xc0000409` crashes to the cylindrical ellipticity/global modal tracking path. The first source audit did not establish a raw `new/delete`, double-free, or obvious RAII ownership defect. Version 1.2.91 therefore applies only conservative boundary hardening rather than speculative algorithm changes:

- IC(0) explicitly invalidates its state before every factorization attempt, preventing stale `valid=true` state after a failed rebuild.
- Cylindrical modal-field validation now checks size arithmetic before multiplication and validates element connectivity again at evaluation boundaries.
- The modal-field locator checks grid-size arithmetic and its integer element-index range.
- Modal-state shape construction and cross-mesh FEM vector checks use overflow-safe size arithmetic.
- Cylindrical global assembly checks local-DOF multiplication and avoids overflowing reserve-size calculations.

No change was made to the eigenvalue algorithm, modal ordering, ellipticity equations, or physical model. These checks are intended to convert silent memory-corruption risks into deterministic exceptions/validation failures. The two original `0xc0000409` failures remain a separate unresolved reproduction until the Windows failing process is localized with a debugger or a sufficiently fast sanitizer/debug build.

## 1.2.93 numerical portability hardening

The cylindrical modal tracking audit exposed a numerically unstable 2x2 principal-angle calculation: the smaller singular value was recovered from `trace^2 - 4 det^2`, which loses the small eigenvalue gap when both singular values are close to one. This could make a self-comparison of an otherwise identical two-dimensional modal subspace fall below the `1 - 1e-9` contract on one compiler/runtime while passing on another. The implementation now forms the eigenvalue gap with `hypot` on the symmetric matrix `C^T C` and recovers the smaller singular value from `|det(C)| / sigma_max`.

The regression `modal_subspace_similarity_contract` exercises a nearly coincident principal-angle case and requires the measured overlap to retain the intended 1e-8 separation. The cylindrical angular-quadrature audit was also made independent of eigen-solver roundoff for its primary convergence gate: it compares the assembled stiffness and mass operators directly, while the modal-frequency comparison remains a bounded secondary check. The targeted ellipticity solve uses a 300-iteration solver budget, which is sufficient for the audited case and avoids spending thousands of iterations when a test contract has already converged.

## 1.2.95 modal-subspace conditioning hardening

The two-dimensional modal-subspace path now uses two-pass modified Gram-Schmidt in the mass inner product. The second projection pass is deliberate: nearly degenerate modal families can lose orthogonality during the first pass even when the exact residual is well defined. No arbitrary residual cutoff was introduced, because a small residual can represent a physically valid but closely aligned second direction. Non-finite intermediate inner products and normalization factors now fail the local orthogonalization step deterministically.

The 2x2 principal-overlap result is explicitly clamped to the mathematical interval `[0,1]` after the stable singular-value reconstruction; this only removes floating-point excursions above one. `modal_subspace_conditioning_contract` covers an almost linearly dependent pair of bases spanning the same plane and requires unit principal overlap. The existing `modal_subspace_similarity_contract` remains the near-unit principal-angle accuracy contract.

## 1.2.95 sparse eigensolver audit

The sparse generalized eigensolver now normalizes K and M independently by their maximum diagonal scales before inverse/subspace iteration. The generalized problem is invariant under these positive coefficient scalings; the physical eigenvalue is restored by the K/M scale ratio. This removes a previously untested dependence of PCG and the dense Ritz problem on absolute coefficient magnitude. A new `sparse_solver_scaling_contract` exercises scale factors `1e-12`, `1` and `1e12`.

The internal M-orthogonalization used by the sparse solver now uses two-pass modified Gram-Schmidt and keeps the same normalized mass metric for both projection and norm evaluation. The small Jacobi Ritz solver uses a relative, matrix-scale-aware off-diagonal convergence criterion instead of a fixed absolute `1e-13` threshold. The generalized residual is evaluated in the normalized problem, which is mathematically equivalent to the physical normwise backward error.

Targeted 1.2.95 solver contracts pass, including the ordinary sparse eigensolver contract, production contract, RCM initial-vector contract, the new scaling contract and both modal-subspace contracts. A complete build/CTest run was not declared green in this environment because the full build exceeded the available execution window.

## 1.2.96 sparse Ritz/PCG robustness hardening

The sparse eigensolver now treats the post-Jacobi Ritz vectors as a checked numerical result rather than assuming that an M-orthonormal input basis remains perfectly M-orthonormal after the dense rotation. Each returned vector pair is checked with the normalized mass matrix; convergence requires the maximum mass-orthogonality error to remain below `1e-8`. The result also exposes the maximum normwise eigenpair residual separately from the legacy aggregate field.

The residual check no longer skips near-zero eigenvalues. The normwise backward-error denominator remains well defined for `lambda = 0`, so numerical zero/near-null modes are now required to satisfy the same residual contract instead of being implicitly accepted by an eigenvalue floor. This is important for nearly singular FEM systems where a small eigenvalue may be physical rather than numerical noise.

PCG tolerance is tied to the requested outer tolerance (`min(1e-10, 0.1*tolerance)`, bounded below by `1e-14`) instead of using an unconditional `1e-8`. This prevents a loose linear-solve stopping rule from becoming the accuracy bottleneck while the Ritz residual is still being used as the final convergence gate.

The IC(0) regularization shift remains a solve-only spectral shift: the PCG operator is `K + delta M`, which preserves generalized eigenvectors exactly and shifts eigenvalues by `delta` only inside the linear solve. The selected shift is now reported for diagnostics.

The sparse matrix symmetry validator was also hardened to check both triangular directions. A one-sided lower-triangular coefficient is now rejected instead of being silently accepted as symmetric. Targeted clustered-spectrum, near-null conditioning, scaling, production, and sparse-mesh robustness contracts pass in the local GCC build.
