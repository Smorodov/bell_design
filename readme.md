# BellOptimizer 1.2.140

BellOptimizer is a dependency-light C++17 application for bell geometry optimization, modal analysis and audio synthesis. It uses a meridional Q4 mesh with a vector-valued cylindrical Fourier FEM for non-axisymmetric and elliptic structural analysis. It is not a full 3-D production FEM.

## Current status

- The Western `western_five_partial` optimizer remains the production optimization path.
- Direct cylindrical vector FEM is available for ellipticity-aware analysis and validation.
- The Chinese-bell two-tone work is a detached research module under `research/chinese_bell/` and is **not part of the normal build**.
- The Chinese research handoff is `research/chinese_bell/docs/CHINESE_RESEARCH_HANDOFF_1.2.126.md`.
- The research module preserves the experimental code, tests, literature and negative results so future work can continue without repeating the documented mistakes.
- Enable it explicitly with `-DBUILD_CHINESE_RESEARCH=ON`; normal builds leave it uncompiled.
- The current reduced Chinese model has not yet established a physically identified `(m,n)a/(m,n)b` family with robust reciprocal A/B strike selectivity.
- CalculiX is an external differential-validation reference, not a production dependency.

## Directional acoustics

The opt-in `bell synth-contact` command remains the external translating-striker diagnostic. The new `bell synth-clapper` command uses direct cylindrical vector FEM (real Fourier harmonics m=0..3), actual inner-wall modal fields, a spherical steel head on a rigid pendulum, conservative Hertz contact and the existing complex outer-wall far-field radiation coefficient. Both paths write peak-normalized mono WAV previews, not calibrated SPL or validated room-acoustic predictions. Modal spectral convergence, contact-curvature calibration and room acoustics remain open. See `docs/CONTACT_DRIVEN_AUDIO_PIPELINE.md`, `docs/EXCITATION_AUDIT.md`, `docs/LOCALIZED_STRIKE_AZIMUTHAL_AUDIT.md`, `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md` and `docs/INTERNAL_CLAPPER_PHYSICAL_MODEL.md`.


## Internal clapper foundation (1.2.139)

Version 1.2.139 added the true inner-wall modal export and spherical-tip closest-point geometry. Version 1.2.140 builds on that foundation with the first coupled physical-pendulum, internal Hertz-contact and WAV synthesis path. See `docs/INTERNAL_CLAPPER_PHYSICAL_MODEL.md`.

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
bell synth-contact [geometry.json] [--output file.wav] [--strike-azimuth-deg deg] [--listener-polar-deg deg] [--listener-azimuth-deg deg]
bell synth-clapper [geometry.json] [--output file.wav] [--tip-radius m] [--tip-length m] [--initial-angular-velocity rad/s]
bell --help
bell --version
```

Use `--workdir <dir>` to select the directory containing geometry, target and synthesis configuration files. See `docs/USER_GUIDE.md` for file formats and complete workflows.

## Source map

- `src/core/` — optimizer, shared domain types and modal tracking.
- `src/models/` — geometry priors and profile models.
- `src/fem/` — cylindrical Q4 elements, global assembly, sparse eigensolver and modal evaluation.
- `src/math/` — numerical helpers and the still-used reduced-order coordinate-morphing fallback.
- `src/excitation/` — production meridional influence, localized strike, external-striker Hertz contact, inner-clapper geometry and coupled pendulum/contact dynamics. Chinese-bell research excitation is detached under `research/chinese_bell/`.
- `src/synthesis/` — modal audio synthesis and acoustic transfer.
- `src/ccx/`, `src/validation/` — external/reference validation only.
- `tests/` — mathematical contracts, regression tests and focused convergence checks.

## Documentation

- `docs/CODEBASE_AUDIT.md` — cleanup decisions and release validation evidence.
- `docs/ARCHITECTURE.md` — authoritative module boundaries and data flow.
- `docs/MATHEMATICAL_AND_ALGORITHMIC_FOUNDATIONS.md` — current equations and numerical methods.
- `docs/VALIDATION.md` — test inventory, evidence and known limitations.
- `docs/DEVELOPMENT_PLAN.md` — ordered next steps and release gates.
- `docs/VALIDATION_THRESHOLDS.md` — normative acceptance criteria, thresholds and evidence states.
- `docs/EXCITATION_AUDIT.md` — excitation model contracts and limitations.
- `docs/CONTACT_DRIVEN_AUDIO_PIPELINE.md` — geometry-to-WAV workflow, units and validation boundaries.
- `docs/CLAPPER_GEOMETRY_AND_KINEMATICS_AUDIT.md` — internal spherical-head clapper design audit.
- `docs/INTERNAL_CLAPPER_PHYSICAL_MODEL.md` — equations, assumptions, CLI parameters and validation limits for `synth-clapper`.
- `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md` — authoritative acoustic-field and stereo design gate.
- `docs/MODAL_STATE_TRACKER.md` — diagnostic continuation tracker.
- `docs/MODAL_FIELD_COMPARISON.md` — full-field cross-mesh and cross-Fourier comparison contract.
- `research/chinese_bell/README.md` — detached Chinese-bell research module.
- `research/chinese_bell/docs/CHINESE_RESEARCH_HANDOFF_1.2.126.md` — handoff and negative-result record.
- `research/chinese_bell/docs/CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md` — consolidated experimental conclusions.
- `research/chinese_bell/docs/CHINESE_BELL_LITERATURE.md` — maintained bibliography.
- `docs/TRADITION_AND_CASTABILITY.md`, `docs/TUNING_PROFILES.md` — geometric priors and reference targets.

The design rule is one authoritative implementation per physical mechanism. Experimental paths remain explicitly isolated until their validation gates pass.

## Internal clapper audio synthesis (1.2.140)

A separate `bell synth-clapper` path now models an internal suspended tongue with a spherical steel head. It uses a planar rigid-pendulum equation, the actual inner-wall Fourier modal fields, equal-and-opposite Hertz contact forces, a conservative impact-energy audit, and the outer-wall directional far-field radiation model to create a mono WAV. The earlier `synth-contact` command remains the external translating-striker diagnostic.

Example:

```text
bell synth-clapper bells/best_geometry.json --output audio/bell_clapper.wav --clapper-mass 5.0 --com-length 0.46 --tip-length 0.52 --tip-radius 0.035 --initial-angular-velocity 3.2
```

The defaults are illustrative rather than measured clapper data. The default pivot is placed 30 mm above the inner crown point, the default inertia is the point-mass estimate `m*l_cm^2`, and the Hertz coefficient assumes a locally flat wall. Supply actual pivot coordinates, inertia, head radius and initial swing conditions for a particular instrument. The first version models one impact through first separation; contact curvature, damping/plasticity, repeated strikes, room acoustics, calibrated SPL and stereo are not yet included. See `docs/INTERNAL_CLAPPER_PHYSICAL_MODEL.md` for equations, assumptions and validation limits.
