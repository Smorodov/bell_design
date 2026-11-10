# Validation and evidence — 1.2.172

## Build contract

C++17, standard library only for the normal application. CalculiX is optional. Configure with `BUILD_SMOKE_TESTS=ON` for test executables and CTest registration. `BUILD_CHINESE_RESEARCH=OFF` is the default.

## Maintained test groups

- **Core/domain:** input safety, material presets, geometry/profile, model-layer separation and castability.
- **Sparse/modal solver:** sparse matrix/eigensolver contracts, eigenpair integrity, solver robustness/scaling, cylindrical Q4 assembly, axis regularity and modal subspace comparison.
- **Excitation/contact:** localized strike, inner-surface modal export, hammer/clapper geometry, pendulum/contact energy, modal release and decay.
- **Acoustics/stereo:** pressure scaling/phase, directional-field rotation, listener geometry, finite-distance radiation and real-FEM receiver integration.
- **End-to-end:** CLI audio smoke tests and selected geometry/mesh convergence checks.
- **Diagnostic continuation:** modal state tracker and cross-mesh field comparison remain test-only and are not part of the optimizer.

The three one-shot 1.2.160–1.2.163 result-generating audit tests were removed after their results were consolidated into maintained analysis documents. Their conclusions are not promoted into physical acceptance criteria.

## 1.2.164 cleanup validation

- All CMake targets built successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`.
- 85 CTest cases are registered after removing three one-shot result-generating audits.
- Focused core/FEM/contact/acoustic group: 24/24 passed.
- End-to-end contact/clapper audio and diagnostic modal-tracking group: 6/6 passed.
- A final focused re-run of the core smoke, finite-distance, clapper and modal continuation cases passed 10/10. This batch overlaps the two groups above; do not add the batch counts as unique tests.
- Full 85-test CTest suite was not run; do not infer full-suite status from the focused groups.
- Windows/MSVC was not tested.

## Evidence semantics

A passing numerical contract establishes only its tested property. It does not imply physical calibration, historical authenticity, correct absolute SPL or agreement with a real bell recording. The finite-distance and far-field kernels remain approximate.

## Chinese-bell status

The module is detached and default-off. Structural solver contracts are numerically tested, but the current reduced model does not establish a physical `(m,n)a/(m,n)b` family or reciprocal A/B strike selectivity. See `research/chinese_bell/docs/CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`.

## Rules

- Do not weaken thresholds to make a candidate pass.
- Do not treat frequency rank as modal identity.
- Do not mix Western and Chinese objectives.
- Do not promote diagnostic acoustic proxies into physical gates.
- Do not claim measured accuracy without an independent reference.


## 1.2.165 source-cleanup validation

- Enabled `-Wall -Wextra -Wpedantic` for the CLI target (MSVC `/W4`); removed the two unused helpers found by this diagnostic pass.
- Full configured build completed successfully; no compiler warnings remained in the captured final build.
- Focused regression batch: 10/10 passed (14.47 s). The batch covers core smoke, finite-distance radiation and real-FEM transfer, modal continuation/subspace, clapper dynamics/evaluation/objective and diagnostic tracker continuation.
- The full 85-test suite was not run in this pass. Windows/MSVC was not tested.


## 1.2.166 — WAV writer and shared radiation validation

- Added `wav_writer_contract`: validates canonical RIFF/WAVE PCM16 headers for mono/stereo, channel interleaving, sample rate, block alignment, clipping and signed endpoints; rejects non-finite samples, mismatched stereo buffers, invalid sample rates and an unwritable destination without creating a partial file.
- All three synthesis workflows now use the same `src/audio/wav_writer.*` serialization implementation.
- Added `directional_radiation_equivalence`: direct azimuthal quadrature vs Fourier/Bessel harmonic field, comparing complex amplitudes over 3 frequencies, 3 polar angles and 17 azimuths per condition. Observed max relative complex difference `2.99e-14`, max absolute difference `2.02e-18`.
- Both production far-field paths (`synth-contact` mono and clapper/stereo) now use the harmonic field builder; the direct integral is retained as a test/reference implementation only. Finite-distance remains a distinct approximate kernel.
- Full configured CMake build completed successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`.
- 87 CTest cases are registered in 1.2.166 (85 from the consolidated baseline plus the two new contracts).
- Focused 10-test batch passed 10/10: WAV writer, complex radiation equivalence, field rotation, listener geometry, finite-distance and real-FEM receiver, hammer contact, clapper dynamics, contact audio smoke and clapper audio smoke.
- Full CTest suite was not run. Windows/MSVC was not tested.


## 1.2.167 — directional-field covariance and reuse

- Added no-copy `evaluateRotated(theta, frameRotation)`; it is algebraically equivalent to a materialized coefficient rotation followed by evaluation.
- Expanded `acoustic_directional_field_rotation_contract` with joint-frame covariance, additive rotation composition, multi-angle no-copy equivalence and base-coefficient immutability.
- Full configured CMake build completed successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`.
- Focused acoustic/stereo batch: 3/3 passed (`acoustic_directional_field_rotation_contract`, `directional_radiation_equivalence`, `stereo_listener_geometry_contract`).
- The full CTest suite was not run. Windows/MSVC was not tested.


## 1.2.168 — joint stereo coordinate covariance

- Added `stereo_coordinate_covariance_contract`. It checks common z-axis rotation of listener centre/ear coordinates, preservation of ear distances/polar angles/interaural geometry/propagation delays, invariant strike-to-listener relative azimuth, and complex left/right transfer invariance using no-copy field-frame evaluation.
- Full configured CMake build completed successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`; `bell --version` reports `1.2.168`.
- The seven focused acoustic/stereo/contact tests passed 7/7.
- All 88 tests are registered. A full-suite run passed tests 1–69, then exceeded the execution window while running `cylindrical_ellipticity_convergence` (#70). A continuation passed tests 71–75 and timed out during `cylindrical_unified_ellipticity_audit` (#76). A separate run passed tests 77–88 (12/12). Thus 86 tests are confirmed passing; tests #70 and #76 remain unverified due to runtime limits, not recorded as failures.
- Windows/MSVC was not tested.


## 1.2.170 — clapper/modal frame covariance and polar-angle audit

- Added `clapper_stereo_frame_covariance_contract`: production internal-clapper contact integration under a common rotation of an exactly degenerate cosine/sine modal pair and strike azimuth; compares contact times, force history, peak force, contact point and release modal state, then verifies complex modal radiation covariance at a rotated listener azimuth.
- Added `stereo_polar_angle_approximation_contract`: compares centre-polar far-field sampling with direct per-ear polar/azimuth evaluation and checks the polar-specific Fourier/Bessel field against direct quadrature.
- On the synthetic 640 Hz multi-harmonic fixture, normalized centre-polar shortcut errors were 61.22%, 29.30%, 14.20% and 6.97% at 0.75, 1.5, 3 and 6 m. The polar-specific builder differed from direct quadrature by less than `2e-11` of the reference peak. This is a fixture-specific diagnostic, not a general physical error bound.
- No production contact, radiation, FEM or optimizer formulas were changed; the centre-polar shortcut remains in production and is now explicitly documented as an approximation requiring a polar-dependent follow-up.

- Full configured CMake build completed successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`; `bell --version` reports `1.2.170`.
- Focused 10-test stereo/acoustics/clapper regression batch passed 10/10.
- 90 CTest cases are registered. All tests except the two long ellipticity audits were confirmed passing (88/88 of the non-long tests); `cylindrical_ellipticity_convergence` and `cylindrical_unified_ellipticity_audit` each exceeded a 180-second execution window without producing a test result. They are unverified, not recorded failures. Windows/MSVC was not tested.


## 1.2.170 — polar-dependent directional field

- Added `polar_directional_field_interpolation_contract`: a 91-node global audit grid, 40 complex comparisons against independently built exact-polar fields, local yaw-envelope interpolation at four listener distances and five yaw values, and azimuthal frame-rotation commutation. The global-grid fixture error was 0.442% of peak; the local-grid yaw-envelope error was 0.365% of peak.
- Production stereo far-field now builds one polar family per mode and evaluates the two ear polar/azimuth directions from it. Mono retains one exact-polar field. FEM/contact/modal time histories are still shared.

- Production stereo smoke test: `synth-clapper bells/best_geometry.json --stereo --duration 0.25 --listener-distance 3.0` generated a valid 48 kHz, 16-bit, two-channel PCM WAV using 32 modes. Contact duration was 1.04 ms; reported conservative energy drift was 3.5e-05. Wall time for the full command was 11.88 s in this Linux Release environment (includes geometry/modal/contact work; not a radiation-only benchmark).
- Full configured CMake build completed successfully with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`; CLI reports `BellOptimizer 1.2.170`.
- CTest registered 91 tests. 88 tests were confirmed passing. `cylindrical_ellipticity_convergence` and `cylindrical_unified_ellipticity_audit` were excluded from the broad run because of their previously observed long runtimes; `cylindrical_modal_field_fem_matrix` exceeded the execution time limit in its separate tail run and remains unverified. Six directly relevant stereo/radiation tests passed 6/6.


## 1.2.171 — real-FEM polar interpolation contract

- Added `real_fem_polar_interpolation_contract`, using the production cylindrical FEM modal evaluator rather than a synthetic modal fixture. It covers 3 cases (17-node baseline, 25-node refinement of the same profile, and a second 17-node bell profile), 54 returned modes spanning 79.55–735.10 Hz, 1/3/5 m receiver distances, five yaw values (-90/-45/0/+45/+90 degrees), both ears, and 1,620 complex comparisons.
- Worst per-mode peak-normalized complex interpolation error: `0.000118175` (0.01182%). Maximum phase error above 10% of each mode's peak: `0.00126688` degrees. Both pass the predeclared limits of 1% and 1 degree.
- 162 independent direct-quadrature spot checks agreed with exact-polar Fourier/Bessel evaluation to relative error below `2e-5`.
- This validates interpolation and the harmonic representation numerically for the tested real-FEM cases. It does not validate the physical far-field radiation model against measurements or establish universal bounds for arbitrary geometries.

- Full configured Release build completed with GCC 14.2 / C++17 (`BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`); CLI reports `1.2.171`. Nine directly relevant stereo/radiation tests passed 9/9. A production clapper stereo smoke run generated a valid 48 kHz, 16-bit, 2-channel PCM WAV (0.25 s, 32 modes). 92 CTest cases are registered; the complete CTest suite was not run for this release.


## 1.2.172 — prepared stereo receiver reuse

- Added `PreparedStereoRadiation` to separate immutable per-mode polar-field preparation from receiver orientation queries. The query operation has no FEM, contact-integrator or modal-release dependency.
- Extended the real-FEM polar contract across all three existing geometry/mesh cases: five yaw values at three centre-azimuth offsets per case, both ears, exact-polar Fourier/Bessel comparisons, receiver range/delay checks, stable mode/build counts, and out-of-envelope rejection.
- The 9-test stereo/radiation regression group passed 9/9 after the change.
- Full configured Release build completed successfully with GCC 14.2 / C++17 (`BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`).
- Production `synth-clapper --stereo` smoke generated a 48 kHz, 16-bit, 2-channel PCM WAV, 0.25 s, 32 modes, at 15-degree head yaw.
- 92 CTest cases are registered; the complete suite was not run. Windows/MSVC was not tested.
- This validates reuse of the prepared receiver-query representation. It does not add an interactive UI or persistent cache between CLI processes, and does not validate physical acoustics against measurements.

### 1.2.174 reusable modal-profile bundle

`bell_modal_profile_contract` verifies boundary modal scalar/array round-trip, negative Fourier harmonic labels, identity rejection and reusable receiver preparation after loading the paired profile. `real_fem_polar_interpolation_contract` additionally serializes and reloads four real FEM modes from the 17-node reference geometry, checks exact preservation of inner/outer Fourier boundary traces, and confirms that the corresponding `.acfield` cache can prepare a stereo receiver without surface integration. The structural profile is boundary-complete for the current clapper/radiation path, but it is not a complete FEM restart.


### 1.2.176 standalone synthesis SDK

- `bell_synthesis` compiles as a separate target with no FEM/optimizer/CLI source dependency.
- Installed package exposes `BellSynthesis::bell_synthesis`; external consumer project uses `find_package(BellSynthesis CONFIG REQUIRED)` and links only that target.
- `synthesis_sdk_link_contract` checks profile-loader and clapper-contact symbols through the standalone target.
- The installed static archive was inspected for optimizer/FEM namespace symbols; none were found.
- This verifies package/link isolation, not a full high-level runtime audio-render workflow.


## 1.2.184 Stage 6 focused verification — 2026-11-10

- Fresh CMake configure succeeded with GCC 14.2, Linux, Release, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`.
- `bell_core`, `bell_synthesis`, and the new `round_solver_full_volume_probe` target built successfully. The `round_optimizer_integration_probe` target also built and ran successfully (`modes=26`, `max_abs_projection=0.0660095`, `modal_family_groups=15`).
- Focused regression set passed 5/5: `cylindrical_modal_family_tracking`, `modal_subspace_similarity_contract`, `modal_subspace_conditioning_contract`, `clapper_modal_spectrum_family_audit`, and `modal_family_adapter_contract`.
- Full configured all-target build was attempted twice but exceeded the available execution window while compiling remaining CLI/test targets; no compiler diagnostic was emitted before termination. Therefore a complete all-target build is **not** claimed. Windows/MSYS2 UCRT64 was not tested in this run.
- Stage 6 cross-solver research probe: 61/61 geometries passed at 15 nodes (305/305 first-five matches); first ten geometries passed at 23 and 31 nodes (50/50 each). These use diagonal row-sum mass weights and are not a consistent-mass or physical-acoustic validation.
