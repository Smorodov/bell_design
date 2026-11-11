# BellOptimizer — project guide

This is the single reference for supported workflows, project boundaries, configuration and synthesis. Current validation evidence and known gaps are maintained separately in [VALIDATION.md](VALIDATION.md). Development/release procedure is in [DEVELOPMENT.md](DEVELOPMENT.md); Russian command examples are in [RUN_OPTIMIZATION_EXAMPLES_RU.md](RUN_OPTIMIZATION_EXAMPLES_RU.md).

## 1. Scope and current status

BellOptimizer is a C++17 application for circular-bell geometry optimization, cylindrical/Fourier modal analysis, contact simulation and modal audio synthesis. The structural solver is a meridional Q4 cylindrical finite-element model with Fourier harmonics; it is not a general-purpose 3-D FEM package.

The active production line is **round bells with round profiles**. Chinese-bell work remains detached under `research/chinese_bell/`, disabled by default and outside the production optimizer. The development order is round-bell optimizer, then generator; SoundFont/SF3 work is deferred.

The established `frequency_profile` objective remains the default. `clapper_acoustic` is an experimental round-bell objective using direct cylindrical FEM and a demonstrative pendulum/Hertz clapper model. Its parameters are not calibrated to a measured installation.

## 2. Build and test

Windows/MSYS2 UCRT64:

```bat
build.cmd
test.cmd
```

Portable CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SMOKE_TESTS=ON -DBUILD_CHINESE_RESEARCH=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

`build.cmd` builds; it does not run tests. CalculiX is optional and used only for external reference checks. See `DEVELOPMENT.md` for the validation ladder and release checklist.

## 3. Main commands

```text
bell optimize [--solver=ccx|internal] [--workdir DIR] [--evolution FILE] [--targets FILE]
bell synth [file.dat] [--workdir DIR]
bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav] [--workdir DIR]
bell synth-contact [geometry.json] [--output file.wav] [--workdir DIR]
bell synth-clapper [geometry.json] [--output file.wav] [--stereo] [--workdir DIR]
```

`--workdir` changes the working directory before resolving relative paths. Run independent optimization cases in separate work directories: each writes `best_geometry.json` and `best_bell.dat` and may reuse an existing geometry as its initial guess.

Ready-to-run Windows launchers include `run_western_220.cmd`, `run_western_440.cmd`, `run_western_880.cmd`, and the corresponding `run_orthodox_*.cmd`. `run_all_demo_cases.cmd` runs all six. Each case writes to its own directory under `runs/`.

## 4. Optimization and target profiles

Target profiles use JSON schema v2 and describe acoustic targets independently of the geometry prior and optimization strategy. Bundled targets are in `bells/targets*.json`; search settings and `seed_family` are in `bells/evolution*.json`. `profile_id` and display names are descriptive metadata, not scoring switches. Legacy schema v1 is rejected.

The Western minor-third reference is `Hum : Prime : Tierce : Quint : Nominal = 0.5 : 1 : 1.2 : 1.5 : 2`; the major-third variant uses 1.25 for Tierce. These are reference profiles, not universal prescriptions. The Orthodox/Sysoi ratios (`0.51 : 1.26 : 1.64 : 2.19` relative to F2) are an experimental target based on one historical exemplar, not a universal Russian bell law. Sorted-index mode selectors remain provisional proxies; physical family identity is not yet established. Chinese two-tone profiles are not production-supported.

The two-stage round-bell optimizer evaluates the population with a fast internal solver and optionally refines a shortlist with direct cylindrical FEM. The same target profile is used at both stages; fast and accurate scores must not be mixed in elite ranking. Stage-2 `m=0` crossmatching supports the tested axisymmetric correspondence, but shortlist recall over the full search domain is not yet established. Keep this feature opt-in.

Example experimental settings:

```json
{
  "enable_round_two_stage_optimization": true,
  "round_two_stage_shortlist_size": 4
}
```

Only `sorted_index` selectors are currently implemented. Under `frequency_profile` they address sorted modes; under `clapper_acoustic` they address frequency-sorted groups after near-degenerate pair merging. Neither is a calibrated modal-family classifier.

## 5. Geometry, material and physical models

Material is selected by preset name in the evolution configuration. `bell_bronze_reference` is an engineering reference (`E=11000 MPa`, `nu=0.34`, `rho=8.60e-9 tonne/mm^3`), not a certified historical composition. `bronze_aluminum_reference` is a mechanical reference (`E=110000 MPa`, `nu=0.30`, `rho=7.40e-9 tonne/mm^3`). Unknown preset names are configuration errors. Chemical composition is not modeled.

The `frequency_profile` excitation influence is a simplified proxy. The `clapper_acoustic` path uses a planar physical pendulum, spherical Hertz contact and modal coupling. Its similarity scaling assumes lengths and tip radius proportional to `s`, mass `s^3`, pivot inertia `s^5`, gravity-similar angular velocity `s^-1/2`, Hertz stiffness `s^1/2`, and contact time scales `s^1/2`. Initial pendulum kinetic energy therefore scales as `s^4` for zero initial angle. This is a similarity demonstrator, not a universal scaling law for real installations.

At separation, the clapper model records initial energy, residual clapper energy, total energy, bell modal energy and per-mode energy. The transfer fraction is `bell modal energy / initial clapper energy`. The conservation check covers the represented undamped coupled model; it does not include contact damping, plasticity, radiation during contact or unmodeled modes.

## 6. Modal selection, acoustics and WAV level

A fixed number of modes does not represent a fixed absolute-frequency band after bell scaling. `synth-contact` and `synth-clapper` now expand their direct FEM modal request geometrically (24, 48, 96, ... by default) until the requested absolute-frequency ceiling is covered, the configurable modal cap is reached, or a larger solve fails. The last successful lower-band solution is retained if expansion fails. The ceiling is `min(20 kHz, 0.49 * sample_rate)`. Use `--modal-max-modes` to set the cap (24..512); metadata records solve count, final request, highest returned frequency and whether the band was covered. A cap-hit report is not evidence of complete coverage.

The physical modal selector applies absolute-frequency and numerical-integrity checks. It does **not** yet implement calibrated shape/family purity. Eigensolver residual, orthogonality, or acoustic amplitude must not be treated as a substitute for physical modal purity. Near-degenerate modes require subspace-aware treatment.

Both contact-driven WAV paths default to peak-normalized PCM16 at 0.95 for backward compatibility. Use `--normalization none` to write unclipped IEEE float32 WAV with the raw pressure-proxy amplitude retained; metadata records raw peak/RMS proxy and applied gain. This is useful for controlled scale comparisons, but it is still not calibrated pascals or SPL. Peak-normalized WAV peak/RMS do not preserve absolute transferred energy. See `research/scale_audio_audit/README.md` for the four-scale procedure. The radiation result is a linear proxy, not a calibrated microphone prediction; room acoustics, head/pinna filtering and measured receiver calibration are not implied.

## 7. Synthesis SDK

`BellSynthesis::bell_synthesis` is a standalone C++17 library for host applications consuming precomputed modal/acoustic boundary bundles. It excludes optimizer, mesh generation, FEM and CalculiX. Install and consume it with:

```sh
cmake --install build --prefix <prefix>
```

```cmake
find_package(BellSynthesis CONFIG REQUIRED)
target_link_libraries(my_game PRIVATE BellSynthesis::bell_synthesis)
```

The `.bmodal`/`.acfield` data are boundary/modal caches, not full FEM restarts: they do not contain interior eigenvectors or global stiffness/mass matrices. The current SDK does not provide a high-level `renderStrike()` facade, real-time audio callback, streaming sink, C ABI or Unity/C# binding.

## 8. Configuration and output files

- `evolution*.json`: search parameters, objective selection, material preset and optional clapper settings.
- `targets*.json`: schema-v2 acoustic target profile.
- `synthesis.json`: sample rate, oversampling, duration, losses, normalization and optional reverb.
- `best_geometry.json`: optimized geometry.
- `best_bell.dat`: optional solver-produced modal data; not guaranteed with the internal solver.
- `*_metadata.txt` and `*_modes.csv`: synthesis diagnostics, including energy and per-mode values where available.

`target_peak` controls digital presentation level, not physical SPL. Reverb is a post-effect and disabled by default unless configured. The geometry synthesis path and the two contact-driven paths have different physical assumptions; do not compare their output levels as if they shared a calibrated source model.

## 9. Geometry and casting boundaries

Geometry validation checks positive wall thickness, plausible thickness ratios, sampled wall-thickness gradients, valid inner surfaces and Q4 element/Jacobian quality. Mold filling, gating, feeding, solidification and thermal casting are outside the current model.
