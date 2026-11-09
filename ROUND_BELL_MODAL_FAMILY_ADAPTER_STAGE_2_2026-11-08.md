# Round-bell modal-family adapter — stage 2 (2026-11-08)

## Scope

This stage continues from `BellOptimizer-1.2.184-round-bell-optimizer-integration-stage1-2026-11-08.zip`.
It adds a research-only adapter between direct cylindrical vector-FEM modes and the optimizer's future modal-family interface. It does **not** replace the production solver or change the objective.

## Implemented

- `research/round_bell_audit/modal_family_adapter.hpp/.cpp`
  - Reads the full direct-FEM displacement vector in `[node][component][Fourier basis]` layout.
  - Computes per-mode azimuthal-order energy shares from the exported mass-lumped row weights.
  - Groups adjacent modes whose relative frequency gap is at most `1e-5` by default.
  - Computes a subspace strike-projection norm `sqrt(sum(p_i^2))` for each grouped eigenspace.
  - Reports an ambiguous family (`dominant_azimuthal_order=-1`) rather than forcing a label when no harmonic has at least 50% of the group-average energy.
- `tests/modal_family_adapter_contract.cpp` validates a synthetic circular `m=1` degenerate pair, its projection norm, and invariance under an orthogonal rotation of the pair basis.
- The existing round optimizer integration probe now also emits a family CSV.
- CMake registers the adapter contract and links the adapter into the research probe only.

## Results

The adapter contract passes:

```text
modal_family_adapter_contract: PASS
degenerate_pair_projection_norm=5
orthogonal_rotation_invariance=PASS
separate_family_count=2
```

The real round-bell probe passes with the baseline geometry and strike configuration:

```text
round_optimizer_integration_probe: PASS
modes=26
max_abs_projection=0.0660095
modal_family_groups=15
```

The family report is `audit_results/round_optimizer_integration_probe_stage2.csv.families.csv`.
The 26 returned modes formed 15 frequency groups; eleven groups contain two near-degenerate modes and four groups are singletons. Harmonic participation is cleanly concentrated in one azimuthal order for this test geometry and retained basis (`m=0..3`). This is diagnostic evidence only, not a mapping to Hum/Prime/Tertia/Quinta/Nominal.

## Interpretation and limits

1. The combined projection norm is invariant to sign changes and orthogonal basis rotations inside an exactly degenerate, mass-orthonormal eigenspace. The contract test exercises a nontrivial rotation.
2. The family energy uses the exported row-sum/lumped mass weights. It is a practical diagnostic and is not identical to evaluating the full consistent mass quadratic form per harmonic. That distinction must be assessed before using the values as optimizer fitness terms.
3. Frequency clustering is a numerical policy. A fixed relative tolerance may merge distinct modes or split a pair under mesh changes; the next validation should sweep this tolerance and compare modal continuation across meshes/geometries.
4. The dominant azimuthal order is not a western partial identity. Assigning five target partials requires tracking family branches across geometry and mesh, and comparing physically justified frequency/ radiation/strike measures. No partial labels are assigned here.
5. The production objective `western_five_partial`, `src/core/bell_optimizer.cpp`, and production runtime FEM/excitation paths remain unchanged. The CMake changes add research/test targets only. No functional version bump.

## Reproduction

```bash
cmake -S . -B build_modal_family_stage2 -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON
cmake --build build_modal_family_stage2 --target round_optimizer_integration_probe modal_family_adapter_contract -j2
./build_modal_family_stage2/modal_family_adapter_contract
./build_modal_family_stage2/round_optimizer_integration_probe audit_results/round_optimizer_integration_probe_stage2.csv
```

See the adjacent `audit_results_stage2_*.log` files and SHA-256 manifest for the captured run.
