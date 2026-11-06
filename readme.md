# BellOptimizer 1.2.177

BellOptimizer is a C++17 application for bell geometry optimization, modal analysis and audio synthesis. The structural engine uses a meridional Q4 mesh and reduced cylindrical/Fourier formulations; it is **not** a full 3-D production FEM package.

## Current implementation status

- **Production optimizer:** the Western `western_five_partial` objective remains the active optimization path.
- **Structural analysis:** the established axisymmetric solver is used by the production optimizer. Direct cylindrical vector FEM is available for ellipticity-aware evaluation, excitation and validation workflows; it is not a general-purpose 3-D solver.
- **Internal clapper:** `synth-clapper` simulates a rigid pendulum with a spherical steel tip, inner-wall modal contact and post-separation modal decay. This is a parametric listening model, not yet calibrated against a measured bell/clapper pair.
## Standalone synthesis library

The installable C++17 target `BellSynthesis::bell_synthesis` provides profile loading, cached modal-boundary contact, directional radiation and audio primitives without linking the optimizer, FEM solver, mesh generator or CalculiX adapter. It is intended for games and other host applications that consume precomputed `.bmodal` + `.acfield` profiles. See [`docs/SYNTHESIS_SDK.md`](docs/SYNTHESIS_SDK.md) and the external consumer example in `examples/synthesis_consumer/`. This release provides the modular library target and low-level API building blocks; a one-call real-time rendering facade and C ABI are future work.

- **Stereo:** `synth-clapper --stereo` shares FEM, contact and modal time history between channels. Mono far-field uses one Fourier/Bessel directional field; stereo far-field builds one local polar family per mode and evaluates both ears from that shared family. The far-field stereo path precomputes only the local polar interval needed to cover the supported head-yaw range (-90 to +90 degrees), with a maximum 2-degree node spacing. Head yaw then changes only receiver queries inside this immutable family; it does not rebuild the radiation field. Synthetic interpolation checks measured 0.442% peak-normalized error on the global reference grid and 0.365% across the production yaw envelope (five yaw values, four distances). A real-FEM contract in 1.2.171 additionally covers three geometry/mesh cases and 54 modes: worst per-mode peak-normalized complex interpolation error was 0.01182%, with maximum phase error 0.00127 degrees above 10% of each mode peak. These are numerical interpolation results, not a universal physical error guarantee or measured acoustic validation. `PreparedStereoRadiation` separates polar-family construction from repeated listener-yaw queries. The explicit `AcousticFieldCache` library API serializes full-polar, per-mode complex radiation fields and can prepare receiver queries across listener distances/directions without new radiation surface integrals. Cache identity is bound to geometry, material preset and FEM/modal solver fingerprints. Strike/clapper excitation, contact history, modal release amplitudes/velocities and T60 are deliberately excluded. The paired `BellModalProfileBundle` API now adds `.bmodal` boundary-modal data (inner/outer profiles, frequencies, modal masses and Fourier boundary traces) beside the `.acfield` radiation cache. It reconstructs the boundary-only `ModalMode` data required by the current clapper-contact and radiation stages; it is not a full FEM restart and omits stiffness/mass matrices and interior eigenvectors. The CLI does not yet expose cache-in/out flags or an automatic persistent cache manager. A shared modal-pressure kernel applies each ear propagation delay once as a signed relative time shift; it is covered by the stereo phase/ITD contract. `finite-distance` remains an opt-in approximate comparison path. Neither is a calibrated exterior-acoustics solution.
- **Audio output:** all WAV paths use the shared `src/audio/wav_writer.*` PCM16 RIFF/WAVE serializer with mono/stereo layout and input validation.
- **Optimizer integration:** the clapper acoustic objective is implemented as a tested primitive but is **not** wired into production candidate evaluation. No objective switch has been made.
- **Chinese two-tone research:** isolated under `research/chinese_bell/`, disabled in normal builds, and currently `RESEARCH_BLOCKED` because physical modal-family identity and reciprocal A/B strike selectivity have not been established.
- **CalculiX:** optional external differential reference; not a production dependency.

## Build

Windows/MSYS2 UCRT64:

```bat
build.cmd
```

Portable CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Options:

- `BUILD_SMOKE_TESTS=OFF` — build the application and core library without validation executables.
- `BUILD_CHINESE_RESEARCH=ON` — explicitly build the detached Chinese-bell research module (default OFF).
- CalculiX is optional and only needed for the external reference checks.

## CLI

```text
bell optimize [--solver=ccx|internal]
bell synth [file.dat]
bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav]
bell synth-contact [geometry.json] [--output file.wav] [--strike-azimuth-deg deg]
bell synth-clapper [geometry.json] [--output file.wav] [--stereo]
                 [--stereo-radiation far-field|finite-distance]
bell --help
bell --version
```

Use `--workdir <dir>` to select the directory containing geometry, target and synthesis configuration files. See `docs/USER_GUIDE.md` for workflows and file formats.

## Source layout

- `src/core/` — domain types, material presets, geometry generation and optimizer.
- `src/models/` — geometry priors and profile models.
- `src/fem/` — Q4 elements, cylindrical assembly, sparse eigensolver and modal evaluation.
- `src/math/` — coordinate transforms, spline and strain utilities.
- `src/excitation/` — localized excitation, contact, internal clapper dynamics and modal decay.
- `src/acoustics/` — directional radiation integral, harmonic field, reusable acoustic-field and boundary-modal profile caches, listener geometry and finite-distance receiver kernel.
- `src/audio/` — shared PCM16 WAV serialization and modal pressure time-domain synthesis.
- `src/synthesis/`, `src/commands/` — audio synthesis and CLI workflows.
- `src/ccx/`, `src/validation/` — optional external/reference validation.
- `tests/` — retained mathematical contracts, regression tests and reproducible convergence checks.
- `research/chinese_bell/` — isolated non-production research implementation and documentation.

## Authoritative documentation

- `docs/ARCHITECTURE.md` — module boundaries and data flow.
- `docs/MATHEMATICAL_AND_ALGORITHMIC_FOUNDATIONS.md` — governing equations and numerical methods.
- `docs/CLAPPER_MODEL.md` — clapper/contact/audio model, assumptions and limits.
- `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md` — current acoustic and stereo architecture/results.
- `docs/MODAL_ANALYSIS.md` — modal identity, continuation and convergence evidence.
- `docs/VALIDATION.md`, `docs/VALIDATION_THRESHOLDS.md` — test inventory and acceptance policy.
- `docs/DEVELOPMENT_PLAN.md`, `docs/DEVELOPMENT_WORKFLOW.md` — next gates and workflow.
- `docs/USER_GUIDE.md` — user-facing commands and file formats.
- `docs/CODEBASE_AUDIT.md` — consolidation decisions and known architectural debt.
- `CHANGELOG.md` — concise milestone history.
- `research/CONSOLIDATED_NUMERICAL_RESULTS.md` — durable numeric findings replacing generated CSV/WAV bundles.

Generated WAVs, raw sweep CSVs, build directories, binaries and per-version release snapshots are intentionally not included in the source archive. Generate audio locally with the CLI when listening is needed.
