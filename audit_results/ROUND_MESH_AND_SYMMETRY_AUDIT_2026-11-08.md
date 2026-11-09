# BellOptimizer 1.2.184 — round-bell mesh and symmetry audit

Date: 2026-11-08  
Source checkpoint: `BellOptimizer-1.2.184-round-bell-preintegration-audit-2026-11-08.zip`  
Scope: round bells with circular profiles only. Chinese-bell research is excluded.

## Environment

- Linux x86_64
- GCC 14.2.0
- CMake 3.31.6
- C++17, Release
- `BUILD_SMOKE_TESTS=ON`
- `BUILD_CHINESE_RESEARCH=OFF`
- No production implementation files were changed.
- Seven audit-only executables were added under `research/round_bell_audit/` and conditionally registered in CMake when smoke tests are enabled.

## 1. Reproduction of the checkpoint mesh experiment

The original `tests/mesh_convergence.cpp` was built from the checkpoint and run with profile resolutions 15, 21, 31, 41, 51. Its generated CSV matched `audit_results/mesh_convergence.csv` byte-for-byte.

SHA-256 for both CSV files:

`f7b3ebb9a075a4eaa72d10ae84218a9eda0c09a65aecefcb94b11c3c3b5ede58`

Fundamental frequencies (Hz):

- 15: 443.456615334
- 21: 439.518727953
- 31: 437.353028142
- 41: 436.588435932
- 51: 436.234725303

The prior 41→51 frequency change is 1.013361%. This reproduces the old experiment, but it is not by itself a full mesh-independence proof.

## 2. Profile refinement 41/51/61/71/81

The audit harness used the same fixed reference geometry and extended the sequence to 61, 71 and 81 profile nodes. It compares the first ten modes using the existing nodal field samples and reports frequency and shape differences.

Fundamental frequencies (Hz):

- 41: 436.588435932
- 51: 436.234725303
- 61: 436.043009756
- 71: 435.927638667
- 81: 435.852896968

Adjacent-level results, worst over modes 1–10:

| Refinement | Max frequency change | Min absolute shape correlation | Max normalized shape error |
|---|---:|---:|---:|
| 41→51 | 1.013361% | 0.997347 | 7.279704% |
| 51→61 | 0.553990% | 0.998951 | 4.579339% |
| 61→71 | 0.334915% | 0.999483 | 3.216085% |
| 71→81 | 0.217598% | 0.999705 | 2.430472% |

## 3. Wall-thickness refinement

At fixed profile resolution 31, the thickness mesh was refined from 1 to 4 layers. Fundamental frequencies (Hz):

- 1 layer: 437.353028142
- 2 layers: 429.380161789
- 3 layers: 427.623611495
- 4 layers: 426.982337327

Worst over modes 1–10:

- 1→2 layers: frequency change 3.208110%
- 2→3 layers: frequency change 0.965803%
- 3→4 layers: frequency change 0.386861%, minimum shape correlation 0.999694, maximum normalized shape error 2.475358%

Conclusion: a single element through wall thickness is not sufficient for a final mesh-independence claim. The first frequency changes materially as thickness resolution is increased.

### Multi-geometry thickness refinement

A further research sweep covered nine selected geometries (four boundary-oriented profiles and five seeded profiles) at 41 and 81 profile nodes, with 1/2/3/4 wall layers. The rank-based diagnostic showed the expected large 1→2 change (up to 3.57% at 41 nodes and 3.96% at 81 nodes), then settled to at most 0.526% and 0.559% respectively on 3→4 layers.

To avoid confusing mode-rank changes with physical family changes, a second comparison used the full-volume cylindrical mass metric and a one-to-one assignment for the five lowest modes on three edge geometries at both 41 and 81 profile nodes:

| Profile nodes | Transition | Max matched frequency change | Minimum full-field overlap | Common volume |
|---:|---|---:|---:|---:|
| 41 | 1→2 layers | 3.5693% | 0.988792 | 1.0 |
| 41 | 2→3 layers | 0.8983% | 0.997820 | 1.0 |
| 41 | 3→4 layers | 0.3485% | 0.999449 | 1.0 |
| 81 | 1→2 layers | 3.7169% | 0.979393 | 1.0 |
| 81 | 2→3 layers | 0.9905% | 0.991002 | 1.0 |
| 81 | 3→4 layers | 0.3851% | 0.997324 | 1.0 |

This supports thickness convergence for the five lowest modes on the selected edge geometries. It is not yet a 61-geometry thickness sweep.

## 4. Deterministic geometry sample and modal convergence

`research/round_bell_audit/round_geometry_dataset.cpp` creates 64 deterministic candidates using seed `0xB3112026`, including four explicitly selected boundary-oriented profiles. It validates and builds 41- and 81-node profile meshes.

- 64 candidates generated
- 61 successfully built at both resolutions
- 3 mesh-generation failures retained in the CSV (indices 7, 10, 28), all due to `Element quality limits exceeded`; they were not silently discarded
- 61 accepted geometries passed to the modal convergence run
- 10 modes solved at 41/51/61/71/81 nodes for every accepted geometry
- 61/61 modal runs completed without failure

The audit harness performs an exact maximum-weight one-to-one assignment of the ten modes between adjacent meshes, using a C++17-portable bitmask dynamic program. This avoids relying only on equal frequency rank when mode order changes.

Worst results across all 61 geometries and ten modes:

| Refinement | Max frequency change | Min shape correlation | Max normalized shape error |
|---|---:|---:|---:|
| 41→51 | 1.573665% | 0.928834 | 37.0497% |
| 51→61 | 0.861101% | 0.987456 | 15.7892% |
| 61→71 | 0.527205% | 0.995939 | 9.00306% |
| 71→81 | 0.344881% | 0.997185 | 7.49838% |

The frequency-change gate of 3% is met on both final refinement steps for the full accepted sample. Shape correlation also improves strongly on the final two steps. The larger early shape errors are retained as diagnostics and must not be hidden.

## 5. Circular symmetry, harmonics and basis invariance

The following existing contracts were built and passed:

- `cylindrical_axis_regularity_contract`
- `cylindrical_global_assembly_contract`
- `cylindrical_global_modes_contract`
- `cylindrical_bell_geometry_modal_tracking`
- `cylindrical_harmonic_family_continuation`
- `modal_family_tracking`
- `modal_state_tracker_contract`
- `modal_subspace_similarity_contract`
- `cylindrical_strain_contract`
- `cylindrical_modal_field_cross_mesh`
- `cylindrical_modal_field_fem_matrix`
- `modal_tracker_fem_continuation`

Observed on the real round-bell modal-tracking test:

- At q=1, the tracked pair relative split is `6.39742e-11`.
- At q=0.95, tracked pair subspace similarity is `0.999855`.
- At q=0.90, tracked pair subspace similarity is `0.998856`.

The modal subspace contract separately exercised 30 deterministic basis transformations (rotation, permutation, sign changes and rescaling) and 30 controlled principal-angle cases.

## 6. Fast-vs-direct FEM calibration on the 61-geometry sample

A second research harness compared the fast axisymmetric solver with the direct cylindrical vector FEM solver using a larger direct-mode candidate pool.

The initial 61-case run using the existing `trackModalFamilies()` reported all 61 candidates ambiguous, Spearman 0.661978, maximum first-five frequency error 95.7432%, and minimum boundary similarity 0.00639781. This result was not a valid measure of physical disagreement: the tracker sets `count = min(reference_modes, candidate_modes, 10)`, so with 8 fast reference modes and 16 direct modes it only considered the first 8 direct candidates.

The follow-up audit compared 8 fast modes against 16 direct m=0 candidates on the same physical mesh. The fast mode was embedded in the common vector layout as `[ur, 0, uz]`; the direct solver exports its complete `[ur, uphi, uz]` field. Both were compared with the full-volume mass metric on the exact same SI mesh. A rectangular one-to-one assignment was solved independently using (a) full-volume mass overlap and (b) outer-boundary trace overlap.

Results over all 61 accepted geometries at 15 profile nodes:

- Full-volume minimum overlap: `0.999999999999`; mean: `0.99999999999995`.
- Outer-boundary minimum overlap: `0.999999999999`; mean: `0.99999999999993`.
- The two independent assignments agreed on all 61/61 geometries.
- Common physical volume coverage: `1.0` for every comparison.
- Maximum relative frequency error across the first eight matched modes: `8.35992e-13` (fraction).
- Spearman correlation of the matched frequency-shape objective: `1.0`.
- Maximum absolute objective difference: `5.1e-13`.
- Descriptive shortlist recall: top-5% = 4/4; top-10% = 7/7.
- The 61-case full-field CSV was reproduced byte-for-byte on a second run.

A temporary research-only patch was then applied to `trackModalFamilies()`: it preserves up to ten reference modes but searches up to twenty candidate modes using the existing exact rectangular assignment routine. Re-running the original boundary-trace probe with this patch produced 61/61 unambiguous candidates, minimum similarity `0.999999999996`, maximum first-five frequency error `8.29386e-13` (fraction), Spearman `1.0`, and top-5%/top-10% recall of 100%. After the experiment, the production source file was restored byte-for-byte to the checkpoint version and the complete build succeeded again. The proposed patch is saved separately; it is not integrated.

A cross-resolution spot-check on nine selected geometries (four explicit boundary-oriented cases and five representative seeded cases) also passed at 41 and 81 profile nodes. At both resolutions, minimum full-field overlap, minimum boundary overlap and common-volume coverage were all `1.0`; no solver or assignment failures occurred.

**Interpretation:** for the shared axisymmetric `m=0` problem on a common mesh, the fast solver and direct vector FEM produce effectively identical meridional modes once the direct candidate pool is wide enough. The earlier poor calibration result is explained by candidate-pool truncation in the existing tracker, not by a demonstrated mismatch of the underlying m=0 fields. The tracker itself has not been modified; a production fix must be separately reviewed and regression-tested. These results do not establish multi-harmonic or full acoustic-objective equivalence.

### Exploratory full multi-harmonic candidate-pool probe (not a round-bell promotion gate)

An additional exploratory check compared the fast m=0 fields against a direct vector-FEM spectrum assembled with harmonics `{0,1,2,3}`. Because non-axisymmetric branches are interleaved with the m=0 families, the direct candidate pool had to be expanded to 96 modes. Three explicitly selected edge geometries were run individually; all eight fast m=0 modes were matched with full-volume overlap 1.0, common-volume coverage 1.0, and maximum frequency errors below `1.2e-13` (fraction). The matched direct-mode indices extended as high as 68.

A nine-geometry batch with this 96-mode multi-harmonic solve did not complete within the execution window. This remains an exploratory result, not a completed sample-wide validation. It also demonstrates why a small fixed candidate pool is unsafe when comparing a reduced m=0 solver with a full multi-harmonic spectrum. The current round-bell production path uses the fast m=0 solver; this probe is retained for later multi-harmonic work.

## 7. Build and CTest status

- Full CMake build: completed successfully at 100%.
- CTest inventory: 99 registered tests.
- Full suite with `cylindrical_ellipticity_convergence` excluded: 98/98 PASS, 0 failures, 102.00 seconds.
- `cylindrical_unified_ellipticity_audit` passed separately in 34.12 seconds.
- `cylindrical_ellipticity_convergence` remains unresolved: the previous targeted run timed out at 60 seconds; a direct isolated attempt did not finish within the 240-second tool window and was terminated. It is NOT a pass and is excluded from the 98/98 result.

## 8. Reproducibility artifacts

- `audit_results/mesh_convergence_repeat_baseline.log`
- `audit_results/mesh_convergence_41_81.log`
- `audit_results/thickness_mesh_convergence_31_1_4.log`
- `audit_results/round_geometry_dataset.log`
- `audit_results/round_geometry_modal_convergence.log`
- `audit_results/round_geometry_multifidelity_probe.log`
- `audit_results/round_geometry_fullfield_match_16.log`
- `audit_results/round_geometry_fullfield_match_16_repeat.log`
- `audit_results/round_geometry_multifidelity_expanded_pool.log`
- `audit_results/round_geometry_multifidelity_restored.log`
- `audit_results/modal_tracking_checkpoint_before_pool_fix.cpp`
- `audit_results/round_geometry_fullfield_match_41_subset.log`
- `audit_results/round_geometry_fullfield_match_81_subset.log`
- `audit_results/round_geometry_thickness_41.log`
- `audit_results/round_geometry_thickness_81_9.log`
- `audit_results/round_geometry_thickness_mass_41_3.log`
- `audit_results/round_geometry_thickness_mass_81_3.log`
- `audit_results/round_geometry_multiharmonic_one_96.log`
- `audit_results/round_geometry_multiharmonic_case1_96.log`
- `audit_results/round_geometry_multiharmonic_case3_96.log`
- `audit_results/round_mesh_targeted_ctest_final.log`
- `audit_results/round_symmetry_harmonic_ctest.log`
- `audit_results/round_rotation_basis_ctest.log`
- `audit_results/round_full_build_resume.log`
- `audit_results/round_full_build_final.log`
- `audit_results/round_full_ctest_excluding_slow_serial.log` (completed 98/98 run)
- `audit_results/round_full_ctest_final_excluding_slow.log` (later repeat attempt; incomplete during test 83, not counted)
- `audit_results/round_ctest_remaining_serial.log` (partial CTest attempt; individual remaining tests were subsequently checked)
- `research/round_bell_audit/round_geometry_dataset.cpp`
- `research/round_bell_audit/round_geometry_modal_convergence.cpp`
- `research/round_bell_audit/round_geometry_multifidelity_probe.cpp`
- `research/round_bell_audit/round_geometry_fullfield_match.cpp`
- `research/round_bell_audit/round_geometry_thickness_sweep.cpp`
- `research/round_bell_audit/round_geometry_thickness_mass_sweep.cpp`
- `research/round_bell_audit/round_geometry_multiharmonic_match.cpp` (exploratory, incomplete sample-wide run)
- `research/round_bell_audit/modal_tracking_candidate_pool_expansion.patch` (temporary, not integrated)
- `build_mesh_repeat/round_geometry_dataset.csv`
- `build_mesh_repeat/round_geometry_modal_convergence.csv`
- `build_mesh_repeat/round_geometry_multifidelity.csv` (original truncated-pool result)
- `build_mesh_repeat/round_geometry_multifidelity_expanded_pool.csv` (temporary patch result)
- `build_mesh_repeat/round_geometry_multifidelity_restored.csv` (checkpoint behavior reverified after source restoration)
- `build_mesh_repeat/round_geometry_fullfield_match_16.csv` and repeat CSV
- `build_mesh_repeat/round_geometry_fullfield_match_41_subset.csv`
- `build_mesh_repeat/round_geometry_fullfield_match_81_subset.csv`
- `build_mesh_repeat/round_geometry_thickness_41.csv`
- `build_mesh_repeat/round_geometry_thickness_81_9.csv`
- `build_mesh_repeat/round_geometry_thickness_mass_41_3.csv`
- `build_mesh_repeat/round_geometry_thickness_mass_81_3.csv`
- `build_mesh_repeat/round_geometry_multiharmonic_one_out_96.csv`
- `build_mesh_repeat/round_geometry_multiharmonic_case1_out_96.csv`
- `build_mesh_repeat/round_geometry_multiharmonic_case3_out_96.csv`
- `build_mesh_repeat/mesh_convergence_41_81.csv`
- `build_mesh_repeat/thickness_mesh_convergence_31_1_4.csv`

Key SHA-256:

- `round_geometry_dataset.csv`: `8959c4ca91ee46a997c0150ea917ebb409292c6d8e188d2dcbca9b1088d41df7`
- `round_geometry_modal_convergence.csv`: `f8002d43b07144855a15c6b11278db8d374fce6c98ebfbcec892e6e136120f63`

Additional SHA-256:

- `round_geometry_fullfield_match_16.csv`: `c3e624b908623bd7810ff4bc3fc24b3b7f3628ce9711e23f66ea673873b55700`
- `round_geometry_multifidelity_expanded_pool.csv`: `e70eccbe5a126b6617bd2d63fd1c381d87e52c919e6fce8b68fbf1c792bfecc6`
- `modal_tracking_candidate_pool_expansion.patch`: `74ca37e7db27fcc28160dfbe8088fd178f54e5abf7f39cffe54e2db7c1d49385`
- `round_geometry_thickness_mass_41_3.csv`: `4149aa31420bac4538a55375bc2a07a3ae2007313a0a685c97a8254160513dc4`
- `round_geometry_thickness_mass_81_3.csv`: `1fb1439cdf7ba0f0181d62c7de38428a48847adf361082d74d8ec49934633ff7`
- `round_geometry_multiharmonic_one_out_96.csv`: `b1d610a0916d22813c6a43e2f4423fee2593a47fb63e6f0d52627621fe5a319f`
- `round_geometry_multiharmonic_case1_out_96.csv`: `7d2bce4981759a66c75d1120ab0a603702de98a0bd0bc8f384356bcd648693f4`
- `round_geometry_multiharmonic_case3_out_96.csv`: `ddf752d548d054dca91fe7e1f73c236a3beb8fa728f54ea44acdc45af2d958c2`
- repeated 15-node full-field CSV: same hash (byte-for-byte reproducibility)
- `round_geometry_fullfield_match_41_subset.csv`: `60c1326023d99f11f9142f2d3f11be27b269505d9939de2bf3d4317f295d2db8`
- `round_geometry_fullfield_match_81_subset.csv`: `387a5f82e586254b7f0f2f8585a7585448ca2a935b12dbed6de85e04ecda35d2`

## 9. Readiness decision

**Status: NOT READY FOR PRODUCTION INTEGRATION.**

Positive evidence: the round-bell profile-frequency convergence criterion is met on the final two refinements for 61 accepted geometries; circular degeneracy and subspace basis-invariance contracts pass; full build succeeds; 98 of 99 registered CTest tests pass in the completed filtered run.

Remaining gates:

1. Resolve or replace the slow `cylindrical_ellipticity_convergence` test and obtain a complete 99/99 run.
2. Extend the combined profile-and-thickness full-field sweep beyond the three edge geometries already checked at 41 and 81 nodes; the 61-geometry profile sweep still uses one thickness layer.
3. Extend the fast-vs-direct comparison beyond m=0 to the relevant multi-harmonic families and near-degenerate subspaces; the m=0 same-mesh field/mass comparison is now validated on 61 geometries.
4. Correct and regression-test the candidate-pool truncation in `trackModalFamilies()` without prematurely enabling production screening. The m=0 top-5% and top-10% recall are both 100% on this sample, but the full multi-harmonic optimizer shortlist recall remains unvalidated.
5. Verify the build and tests on Windows/MSYS2 UCRT64.
6. Keep the experimental metric out of the production optimizer and leave `western_five_partial` unchanged until all gates pass.

No project version bump was made. This audit does not authorize integration.
