# Round-bell modal-family adapter — Stage 3 (2026-11-08)

## Scope

This stage hardens the research-only modal-family adapter introduced in Stage 2. It does not modify the production optimizer, production modal tracking, or the `western_five_partial` objective. Project version remains 1.2.184. The Chinese-bell branch and SF3 work remain out of scope.

## Changes

1. Added strict common-mesh validation across modal vectors. Modes must use the same harmonic basis and node coordinates (within 1e-12); equal node counts alone are no longer accepted as proof of compatible DOF layout.
2. Made harmonic-order magnitude conversion safe for the minimum signed integer edge case and reject malformed bases whose absolute harmonic order exceeds 4096.
3. Extended `modal_family_adapter_contract` with negative tests for mismatched mesh coordinates and unsorted eigenfrequencies, retaining the orthogonal-basis-rotation invariance check.
4. Preserved the existing frequency grouping rule and research-only status. No western partial labels are inferred.

## Validation

Linux/GCC 14.2, C++17, Release, `BUILD_SMOKE_TESTS=ON`.

- Configure: PASS (`audit_results_stage3_configure.log`).
- Build `modal_family_adapter_contract`: PASS.
- Build `round_optimizer_integration_probe`: PASS.
- Build `circular_localized_strike_projection`: PASS.
- CTest targeted set: 2/2 PASS (`audit_results_stage3_ctest.log`).
  - `circular_localized_strike_projection`: PASS, 5.92 s.
  - `modal_family_adapter_contract`: PASS.
- Round geometry probe: PASS; 26 modes, 15 modal groups; strike projection maximum absolute value 0.0660095 (`audit_results_stage3_probe.log`).

A full project build and full CTest suite were not run in this stage. Windows/UCRT64 validation was not run in this stage.

## Interpretation and remaining limitations

The adapter now rejects modal vectors from different meshes instead of silently comparing their DOFs. The synthetic contract verifies the degenerate-subspace projection norm is invariant under an orthogonal rotation of the two-mode basis. This does not yet prove robustness across a sweep of actual bell geometries, mesh densities, or near-degenerate tolerance choices.

The grouping remains based on adjacent frequency gaps. A chain of individually close frequencies can therefore form a wider group; a future stage should test and, if needed, replace this policy with a bounded cluster criterion. Harmonic energy shares are computed per mode and averaged within each group; they are not yet a full subspace-invariant harmonic decomposition for arbitrary mixtures. No mapping to Hum/Prime/Tertia/Quinta/Nominal is authorized yet.

## Next step

Run a controlled parameter sweep over several round-bell geometries and mesh resolutions, recording modal-family stability versus frequency tolerance and harmonic-share threshold. Only after those results should family-to-partial mapping and opt-in optimizer integration be considered.
