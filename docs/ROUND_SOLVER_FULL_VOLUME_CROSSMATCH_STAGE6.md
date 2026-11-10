# Round-bell cross-solver full-volume cross-match — Stage 6

## Status

Research diagnostic only. No production optimizer code, fitness function, candidate rejection, or modal ordering policy was changed. The project version remains **1.2.184**. Scope is circular bells (`ellipticity=1.0`) and the axisymmetric meridional `m=0` block; Chinese-bell work remains detached.

## Why this stage was needed

The earlier `round_geometry_multifidelity_probe` compared the first eight legacy and direct-FEM modes using unweighted outer-wall traces and a square 8-by-8 assignment. On the 61 valid geometries it marked all 61 assignments ambiguous, with Spearman objective-rank correlation 0.661978 and maximum first-five matched-frequency relative error 0.957432. That result did not establish that the first five physical meridional modes disagree: it also forced the eight legacy modes to compete for all eight direct modes, although the direct vector formulation can return additional `m=0` modal content that the legacy meridional formulation does not represent.

## New comparison method

`research/round_bell_audit/round_solver_full_volume_probe.cpp` performs a one-to-one rectangular assignment:

- Reference rows: the first five modes from the legacy axisymmetric meridional solver.
- Candidate columns: up to the first eight modes from direct cylindrical vector FEM, harmonic block `m=0`.
- Mesh: same generated meridional mesh and node ordering. The probe explicitly verifies node correspondence after converting the legacy mesh coordinates from millimetres to the direct model's SI units.
- Similarity: absolute full-mesh meridional displacement overlap using direct-FEM diagonal row-sum mass weights. The radial and axial DOFs are compared; the azimuthal DOF is excluded because it has no counterpart in the legacy meridional solver. This is a lumped-mass diagnostic, not a consistent-mass inner product.
- Assignment: each of the five legacy modes is assigned to a distinct direct-FEM candidate among the eight available modes. Higher legacy modes do not displace one of the five target rows from the assignment.

## Results

| Mesh resolution | Valid geometries | Matched pairs | Minimum overlap | Maximum relative frequency difference |
|---:|---:|---:|---:|---:|
| 15 nodes | 61 | 305 | 0.999999999999 | 8.29e-13 |
| 23 nodes | 10 | 50 | 1.000000000000 | 3.25e-13 |
| 31 nodes | 10 | 50 | 0.999999999999 | 7.55e-13 |

All tested geometries solved successfully; no mesh-ordering or assignment failures occurred. For the 15-node run, the direct-mode indices selected for legacy modes 0–4 were respectively: mode 0 → index 0 in all 61 cases; mode 1 → index 1 in 55 cases and index 2 in 6; mode 2 → index 3 in all 61; mode 3 → index 4 in 33 and index 5 in 28; mode 4 → index 6 in 44 and index 7 in 17. This variation is consistent with direct candidate ordering changes/crossings and is why the assignment must be based on field similarity rather than fixed frequency rank.

## Interpretation

For the tested circular geometries and matched meshes, the first five legacy meridional modes have essentially identical full-mesh radial/axial fields and frequencies in the direct `m=0` FEM result when matched with the rectangular five-to-eight assignment. The earlier poor calibration result is therefore likely dominated by the previous all-eight square assignment and/or boundary-trace matching, rather than a physical discrepancy in these first five meridional modes. This is a strong diagnostic result, not yet a general solver-equivalence certificate.

## Remaining gates

1. Perform an A/B ablation on the same 61 geometries: square 8-by-8 boundary-trace assignment, rectangular 5-by-8 boundary-trace assignment, and rectangular 5-by-8 full-volume mass-weighted assignment. Record assignment indices and all failures.
2. Replace the diagonal row-sum mass approximation with the consistent FEM mass bilinear form, or document and validate the approximation against it.
3. Extend near-degenerate matches to rank-2 subspace comparisons so arbitrary basis rotation does not cause unstable individual-mode labels.
4. Expand 23- and 31-node checks beyond the initial ten geometries, then test cross-mesh interpolation/overlap.
5. Validate the full intended optimizer geometry domain and report shortlist recall before considering any multi-fidelity candidate screening.
6. Keep `western_five_partial` unchanged and keep all diagnostics opt-in until these gates pass.

## Reproduction

From the project root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=OFF
cmake --build build --target round_solver_full_volume_probe --parallel
build/round_solver_full_volume_probe research/data/round_geometry_dataset.csv research/results/round_solver_full_volume_overlap_15node.csv 15
build/round_solver_full_volume_probe research/data/round_geometry_dataset.csv research/results/round_solver_full_volume_overlap_23node_first10.csv 23 10
build/round_solver_full_volume_probe research/data/round_geometry_dataset.csv research/results/round_solver_full_volume_overlap_31node_first10.csv 31 10
```

The output CSVs preserve per-geometry, per-mode matches. The input geometry dataset is retained under `research/data/` as the minimal reproducibility input.
