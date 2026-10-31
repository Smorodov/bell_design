# BellOptimizer 1.2.95

BellOptimizer is a dependency-light C++17 application for bell geometry optimization, modal analysis and audio synthesis. It uses a meridional Q4 mesh with a vector-valued cylindrical Fourier FEM for non-axisymmetric and elliptic structural analysis. It is not a full 3-D production FEM.

## Current status

- The Western `western_five_partial` optimizer remains the production optimization path.
- Direct cylindrical vector FEM is available for ellipticity-aware analysis and validation.
- The localized Chinese-bell strike and two-tone objective are diagnostic/experimental. They are **not connected to the production optimizer**; the tested two-tone excitation gate has not passed robustly under mesh and Fourier-basis changes.
- The modal state tracker is diagnostic. It uses Kalman frequency prediction and shape/subspace association; it does not modify eigenpairs or maintain a persistent multi-hypothesis bank.
- Full cylindrical-FEM modal fields can be compared across different Q4 meshes and nested Fourier bases using symmetric cross-mesh quadrature in a common physical domain and a density-weighted mass inner product. Version 1.2.92 adds a real-FEM continuation matrix for the previously problematic 11th/12th modal family: mesh refinement, Fourier enrichment, and simultaneous changes all retain high two-dimensional subspace similarity. The utility remains diagnostic in `bell_experimental`; it is not wired into the state tracker or production optimizer.
- CalculiX is an external differential-validation reference, not a production dependency.

## Build and test

Windows/MSYS2 UCRT64:

```bat
build.cmd
```

Manual CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Set `-DBUILD_SMOKE_TESTS=OFF` to build only the application and core library. No third-party C++ libraries are required. CalculiX is optional and is only needed for the external reference check.

## Command line

```text
bell optimize [--solver=ccx|internal]
bell synth [file.dat]
bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav]
bell --help
bell --version
```

Use `--workdir <dir>` to select the directory containing geometry, target and synthesis configuration files. See `docs/USER_GUIDE.md` for file formats and complete workflows.

## Source map

- `src/core/` — optimizer, shared domain types and modal tracking.
- `src/models/` — geometry priors and profile models.
- `src/fem/` — cylindrical Q4 elements, global assembly, sparse eigensolver, direct modal evaluation and cross-mesh full-field comparison.
- `src/math/` — spline/numerical helpers and the still-used reduced-order coordinate-morphing fallback.
- `src/excitation/` — meridional influence, localized strike and Chinese strike model. The isolated Chinese pair objective and modal state tracker are compiled into the test-only `bell_experimental` library. The independent 3-D reference mesh builder is in the test-only `bell_validation` library.
- `src/synthesis/` — modal audio synthesis and acoustic transfer.
- `src/ccx/`, `src/validation/` — external/reference validation only.
- `tests/` — mathematical contracts, regression tests and convergence audits.

## Documentation

- `docs/CODEBASE_AUDIT.md` — cleanup decisions and release validation evidence.
- `docs/ARCHITECTURE.md` — authoritative module boundaries and data flow.
- `docs/MATHEMATICAL_AND_ALGORITHMIC_FOUNDATIONS.md` — current equations and numerical methods.
- `docs/VALIDATION.md` — test inventory, evidence and known limitations.
- `docs/DEVELOPMENT_PLAN.md` — ordered next steps and release gates.
- `docs/EXCITATION_AUDIT.md` — excitation model contracts and limitations.
- `docs/MODAL_STATE_TRACKER.md` — diagnostic continuation tracker.
- `docs/MODAL_FIELD_COMPARISON.md` — full-field cross-mesh and cross-Fourier comparison contract, including the 1.2.92 FEM continuation matrix.
- `docs/CHINESE_BELL_ACOUSTIC_CANON.md` and `docs/CHINESE_BELL_STATUS.md` — Chinese-bell design contract and current readiness.
- `docs/TRADITION_AND_CASTABILITY.md`, `docs/TUNING_PROFILES.md` — geometric priors and reference targets.

The design rule is one authoritative implementation per physical mechanism. Experimental paths remain explicitly isolated until their validation gates pass.
