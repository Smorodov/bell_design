# Round-bell optimizer integration — Stage 1

Date: 2026-11-08  
Base: BellOptimizer 1.2.184  
Scope: circular/round-profile bells only (`ellipticity = 1.0`). Chinese-bell research and SF3/SoundFont are out of scope.

## Purpose

This stage builds and runs a focused integration probe across the existing production components:

1. `BellGeometry` + `BellMaterial` → `make_cylindrical_bell_model()`;
2. direct cylindrical vector FEM → `evaluate_cylindrical_bell_modes()`;
3. exported Fourier boundary fields → `modalImpulseCircular()`;
4. ordered frequency/projection evidence → CSV.

The probe is deliberately research-only. It does not change `EvolutionaryOptimizer::evaluateCandidate()`, `computeError()`, `western_five_partial`, or the default optimizer behavior. No production runtime source file was modified.

## Geometry and strike used

A representative valid round profile based on the project's saved best-geometry proportions, canonicalized to `H=1000 mm`:

- `R_mouth = 540.7 mm`
- `t_top = 34.7 mm`
- `t_bottom = 70.8 mm`
- normalized profile knots: `{0, 0.130778, 0.216866, 0.45905, 1}`
- `ellipticity = 1.0`

FEM settings:

- harmonics `{0, 1, 2, 3}` (both cosine/sine terms for nonzero harmonics)
- angular samples: 64
- requested modes: 18 (solver returned 26 after complete-pair/rigid-mode handling)
- maximum iterations: 1600
- tolerance: `1e-6`

Temporary strike settings (not yet production-calibrated):

- center: 70% of outer-profile axial span (`z = 0.700 m` in this canonical mesh)
- meridional Gaussian sigma: 12 mm
- circumferential Gaussian sigma: 8 mm
- azimuth: 0 rad
- direction: surface normal

## Results

`round_optimizer_integration_probe` built and ran successfully:

- direct cylindrical modal evaluation: PASS
- circular localized strike projection: PASS
- modes returned: 26
- maximum absolute strike projection: `0.0660095`
- CSV: `audit_results/round_optimizer_integration_probe.csv`

The existing `circular_localized_strike_projection` contract test was rebuilt and passed:

- 1/1 passed
- elapsed: 3.24 s

The full Windows CTest log supplied by the user independently reports 99/99 passed. The full Linux CTest suite in the preceding continuation audit also reports 99/99 passed.

## Critical integration finding

The direct cylindrical result is not a drop-in replacement for the current optimizer's five-slot modal array.

The first returned frequencies include near-exact pairs, e.g. approximately:

- 81.0558 Hz (two modes)
- 168.6636 Hz (two modes)
- 212.2731 Hz (two modes)
- 291.0565 Hz (two modes)

The local strike projections within these pairs differ substantially. This is expected for a particular cosine/sine basis orientation, but it means that assigning a target amplitude to one raw eigenvector is not rotationally invariant. A small perturbation or a different eigensolver basis could rotate the pair and change the per-mode projection without changing the physical two-dimensional modal subspace.

Also, `ModalMode::fem_harmonics` currently records the complete retained Fourier basis for each mode; it does not identify the mode's dominant harmonic family. Therefore the CSV harmonic column is the retained basis list, not a classification of each eigenmode.

Consequently, directly replacing the legacy `raw_data` in `evaluateCandidate()` and feeding its sorted individual modes into the existing `computeError()` would silently change the meaning of the `western_five_partial` objective. The current objective assumes stable semantic slots (Hum/Prime/Tertia/Quinta/Nominal), while the new FEM output is a sorted set of vector eigenmodes with degeneracies and multiple circumferential families.

## Integration decision

**Do not enable the cylindrical round FEM path in production optimizer yet.** The next required implementation unit is a modal-family adapter between the cylindrical FEM and the existing objective.

Required work, in order:

1. Export per-mode harmonic participation/energy by Fourier order from the eigenvector, not only the list of available basis labels.
2. Identify near-degenerate cosine/sine pairs and construct basis-invariant family quantities (frequency representative and combined strike projection, e.g. subspace norm).
3. Define the exact mapping from the Western five target partials to the appropriate modal families. Do not infer this mapping from raw frequency rank alone.
4. Validate the family mapping on a geometry batch and against the legacy optimizer path.
5. Only then add an opt-in optimizer setting and compare full objective values. Keep it disabled by default until acceptance criteria pass.

## Source-change boundary

Modified:

- `CMakeLists.txt`: adds one `BUILD_SMOKE_TESTS`-guarded research executable target.
- `research/round_bell_audit/round_optimizer_integration_probe.cpp`: new standalone diagnostic source.

Not modified:

- `src/core/bell_optimizer.cpp` / `.h`
- `src/core/bell_types.h`
- `src/fem/cylindrical_modal_evaluator.cpp`
- `src/excitation/localized_strike.cpp`
- production target/objective data, including `western_five_partial`

No functional version bump. This archive is an integration-stage research checkpoint, not a production-integrated release.

## Reproduction

Configure and build:

```bash
cmake -S . -B build_integration_probe -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=OFF
cmake --build build_integration_probe --target round_optimizer_integration_probe -j2
./build_integration_probe/round_optimizer_integration_probe audit_results/round_optimizer_integration_probe.csv
```

Run the strike contract:

```bash
cmake --build build_integration_probe --target circular_localized_strike_projection -j2
ctest --test-dir build_integration_probe --output-on-failure -R '^circular_localized_strike_projection$'
```

