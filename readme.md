# BellOptimizer 1.2.184 — round-bell development line

BellOptimizer is a C++17 application for circular-bell geometry optimization, cylindrical/Fourier modal analysis, contact simulation and modal audio synthesis. It is a reduced-order meridional Q4 cylindrical FEM model, not a general-purpose 3-D FEM package.

## Current status

- Active target: round bells with round profiles. Chinese-bell research is detached and disabled by default.
- `frequency_profile` remains the default objective; `clapper_acoustic` is experimental and its clapper settings are not calibrated to measured hardware.
- Two-stage fast/FEM optimization is opt-in. Tested `m=0` modal correspondence is strong, but shortlist recall over the full search domain is not yet established.
- Physical modal selection checks absolute frequency and numerical integrity; calibrated shape/family purity is not implemented.
- `synth-contact` and `synth-clapper` adaptively expand modal requests toward a common absolute-frequency ceiling; coverage remains bounded by mesh resolution and the configured mode cap.
- Contact-driven WAV files remain peak-normalized by default for compatibility. Use `--normalization none` for unclipped float32 proxy-level output when comparing scales; it is not calibrated SPL. See the scale-audio audit guide.
- `BellSynthesis::bell_synthesis` is a separate SDK target for host applications using precomputed modal/acoustic boundary bundles.

## Build and test

Windows/MSYS2 UCRT64: `build.cmd`, then `test.cmd`.

Portable CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Documentation

- [Project guide](docs/PROJECT_GUIDE.md) — supported commands, architecture, profiles, physics, synthesis and SDK.
- [Validation status](docs/VALIDATION.md) — measured scale/energy results, test status and open gates.
- [Development workflow](docs/DEVELOPMENT.md) — validation ladder, evidence retention and release checklist.
- [Russian optimization examples](docs/RUN_OPTIMIZATION_EXAMPLES_RU.md).
- [Consolidated historical numerical results](research/CONSOLIDATED_NUMERICAL_RESULTS.md).
- [Bell-bronze damping reference](research/BELL_BRONZE_MODAL_DAMPING.md).
- [Scale-dependent audio audit](research/scale_audio_audit/README.md) — four-scale case preparation and raw-level comparison.

## Source layout

- `src/core/` — geometry, materials, optimizer and candidate selection.
- `src/fem/` — cylindrical Q4 assembly, sparse eigensolver and modal evaluation.
- `src/excitation/` — localized strike, hammer contact and internal clapper dynamics.
- `src/acoustics/`, `src/audio/`, `src/synthesis/` — radiation, modal pressure and synthesis.
- `src/commands/` — CLI workflows.
- `tests/` — contracts and numerical audits.
- `research/round_bell_audit/` — reproducible round-bell research probes.
- `research/scale_audio_audit/` — scale fixtures and raw audio comparison workflow.
- `research/chinese_bell/` — detached, deferred research branch.
