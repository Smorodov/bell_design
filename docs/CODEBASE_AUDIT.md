# Codebase audit and consolidation record — 1.2.170

## Consolidation performed

- Removed generated WAVs, per-mode CSV sidecars, metadata dumps and historical listening-output folders from the source archive.
- Removed raw stereo listener-matrix and acoustic mesh-convergence CSV files after transferring their metrics, setup and limitations into `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md`, `docs/MODAL_ANALYSIS.md` and `research/CONSOLIDATED_NUMERICAL_RESULTS.md`.
- Removed duplicate CCX fixture files (the identical `bell.dat` and unused `best_bell.inp`); the core smoke test now shares the single `best_bell.dat` fixture.
- Removed three one-shot test executables/sources used only to produce the 1.2.160–1.2.163 research tables: `acoustic_mesh_convergence_audit`, `stereo_finite_distance_comparison_audit`, and `stereo_listener_matrix_audit`. Their results are retained as analysis. Regression contracts for field rotation, listener geometry, finite-distance radiation and real-FEM transfer remain.
- Removed append-only per-version release and validation snapshots. Their durable milestones are consolidated into `CHANGELOG.md` and `docs/VALIDATION.md`.
- Consolidated architecture, modal tracking, clapper model and stereo decisions into current authoritative documents. Removed stale duplicate audit documents that contradicted later implementation state.
- Kept the Chinese research module isolated and default-off; it is not linked to production optimization.

## Production boundaries retained

- `bell_core` owns production solver/optimizer/domain/excitation/acoustics code.
- `bell_experimental` is created only when smoke tests are enabled and contains diagnostic modal tracking/full-field comparison helpers.
- `bell_validation` contains the independent reference mesh used by its contract test.
- Chinese-bell research remains under `research/chinese_bell/` and is compiled only with `BUILD_CHINESE_RESEARCH=ON`.
- The clapper acoustic objective is not connected to production optimizer candidate evaluation.

## Source review

All translation units under `src/` are explicitly included by the CMake production or validation/experimental targets. No unreferenced production `.cpp` file was identified in the 1.2.163 source inventory. Therefore no production implementation source was deleted merely to reduce the archive size.

## Test policy

Keep tests that protect solver algebra, element/assembly contracts, modal identity, excitation/contact, acoustic phase/rotation, listener geometry and end-to-end smoke paths. Remove one-shot research harnesses only after recording their numerical findings. A test PASS proves only the tested numerical contract; it does not establish experimental acoustic validity.

## Known architectural debt

- `synthesize_contact.cpp` remains a large CLI orchestration unit and is a candidate for later decomposition, but splitting it during a cleanup-only release would increase regression risk.
- The legacy coordinate-morphing path remains active; it is not interchangeable with direct vector FEM.
- Modal tracking responsibilities are deliberately separate; cross-mesh Q4 field comparison is not wired into the Kalman diagnostic tracker.
- The finite-distance acoustic kernel is approximate and opt-in.

## Archive policy

Do not commit build directories, binaries, generated WAVs, raw sweep tables, transient logs or duplicate per-version validation reports. Keep numerical results, assumptions, thresholds and known limitations in maintained Markdown documents.


## Source-level pass — 1.2.165

- Enabled strict compiler diagnostics on the `bell` CLI target as well as the libraries. This exposed two unused internal functions that were not diagnosed while the CLI target had no warning flags: `basisForLabel` in `synthesize_contact.cpp` and `smoothstep` in `bell_synthesis.cpp`. Both were removed after confirming there were no call sites.
- The complete configured CMake build then succeeded with no compiler warnings in the captured build output.
- A focused 10-test regression batch passed: core smoke, finite-distance kernel and real-FEM receiver, modal family/subspace, clapper dynamics/evaluator/objective, and modal tracker/continuation.
- A source scan found no manual `new`/`delete`, `malloc`/`free` use in the production source set. The optimizer owns its solver and seed provider through `std::unique_ptr`; mesh, modal, contact and audio working data are predominantly value/standard-container based. This is a positive ownership signal, not a proof that every lifetime or peak-memory behavior is optimal.
- In 1.2.166, WAV serialization was extracted to `src/audio/wav_writer.*`; direct and harmonic far-field coefficient calculations were extracted to `src/acoustics/directional_radiation.*`. Production mono and stereo far-field paths now share the harmonic builder. The direct quadrature remains a test/reference oracle after a complex-amplitude equivalence contract passed. `synthesize_contact.cpp` still combines CLI parsing and contact/render orchestration and remains a candidate for a later behavior-preserving split.
- `bell_optimizer.cpp` still contains optimizer orchestration and persistence logic in one translation unit. It is not split in this pass because the candidate evaluation and search flow are tightly coupled and no behavior-preserving seam has yet been established.

### Next source refactor gates

1. Split `synthesize_contact` command parsing/orchestration from contact/render services without moving or changing physical equations.
2. Audit optimizer persistence and candidate evaluation boundaries separately; preserve deterministic-seed and objective contracts.
3. Continue stereo coordinate-covariance, immutable-field reuse and phase/timing validation. The current direct-vs-harmonic comparison establishes only the tested far-field integral equivalence.


## Source refactor — 1.2.166

- Added `src/audio/wav_writer.*` as the single PCM16 RIFF/WAVE serializer. It prevalidates sample buffers, checks RIFF size bounds and sample-rate byte-rate overflow, clips normalized samples, writes little-endian headers/data and checks flush/close status. Both mono and stereo layouts are supported.
- Routed `BellSynthesizer::save_wav`, `synth-contact` and `synth-clapper` output through that module; removed the legacy class-local WAV header implementation and command-local duplicate writers.
- Added `src/acoustics/directional_radiation.*`. The direct azimuthal quadrature is retained as a reference evaluator; production mono and stereo far-field paths now use the shared Fourier/Bessel directional-field builder.
- Complex equivalence test passed over 3 frequencies × 3 polar angles × 17 azimuths, with 8192 direct azimuth samples per comparison; maximum relative complex difference `2.99e-14`.
- Added two maintained CTest contracts. The configured CMake build completed and the focused 10-test acoustic/audio/contact regression group passed. The full suite remains unrun.


## Directional-field reuse — 1.2.167

- Added `AcousticDirectionalField::evaluateRotated` to query a rotated coordinate frame directly from the original coefficient set, avoiding coefficient-vector copies for repeated receiver-angle evaluations.
- Extended the rotation contract to cover joint-frame covariance, additive composition, equivalence to materialized rotation and immutability of the base field.
- This is an algebraic capability and test contract. The production stereo command already builds one far-field representation per mode and evaluates both ears from it; this change does not claim to cache fields across separate CLI invocations or across different polar angles.


## Stereo covariance — 1.2.170

- Added a geometry/representation-level common-frame rotation contract covering source-frame rotation, strike-relative azimuth, idealized ear coordinates, complex receiver transfer and propagation delay.
- This is not an end-to-end clapper modal-phase covariance test. The production command still has no global source-orientation abstraction; a full rotated strike/contact/modal integration remains a separate gate.
- The far-field representation remains azimuth-only at the centre polar angle. The ear polar-angle difference is an acknowledged approximation; do not call this a fully 3D directional field.


1.2.170 adds two focused stereo contracts: production clapper/modal frame covariance for an exactly degenerate harmonic pair, and a measured audit of the listener-centre polar-angle shortcut. No production physical equations changed.


## Polar-dependent stereo renderer — 1.2.170

The production far-field stereo branch now builds a local `AcousticPolarDirectionalField` for the full supported head-yaw polar envelope (-90 to +90 degrees), with no more than 2 degrees between nodes. It evaluates both ears from the same per-mode family. The global 91-node builder is retained for validation. Mono and finite-distance paths are unchanged. This changes the acoustic receiver-transfer implementation only; FEM, clapper contact, modal release/time evolution and optimizer objective remain unchanged.


## Real-FEM polar interpolation — 1.2.171

`real_fem_polar_interpolation_contract` exercises the actual cylindrical FEM modal evaluator and the production polar-family builder. Three cases cover 17-node baseline, 25-node refinement, and a second 17-node profile. Across 54 modes (79.55–735.10 Hz) and 1,620 ear/yaw/distance comparisons, the worst per-mode peak-normalized complex interpolation error was 0.01182%; maximum phase error above 10% of each mode peak was 0.00127 degrees. 162 direct-quadrature spot checks passed the `2e-5` relative-error gate. The 2-degree linear complex interpolation is retained. This contract does not prove measured acoustic accuracy or cache/reuse behavior across independent listener updates; the latter is the next audit item.


## 1.2.172 — acoustic receiver-query separation

The far-field stereo receiver transfer is now encapsulated by `Acoustics::PreparedStereoRadiation`. Preparation builds immutable per-mode polar families from modal surface fields; repeated listener-yaw/centre-azimuth queries calculate ear geometry and complex transfers without calling FEM, clapper contact, modal release or the polar surface-integral builder. The object exposes mode count and polar-family build count for deterministic validation. The production clapper stereo command uses this object after contact separation. This is an API-level lifecycle boundary, not a process-persistent cache or interactive UI. The receiver object owns its precomputed fields by value and does not retain references to the caller's modal/profile buffers.


## 1.2.174 reusable modal profile

`src/acoustics/bell_modal_profile.*` owns the versioned boundary-modal profile format. It stores the inner/outer geometry and modal Fourier boundary traces required by the current clapper-contact and directional-radiation modules, with caller-supplied geometry/material/solver fingerprints and an integrity checksum. `BellModalProfileBundle` pairs that structural profile (`.bmodal`) with the existing full-polar radiation cache (`.acfield`). It deliberately excludes strike-specific force, clapper trajectory, contact history, modal release state and PCM. It is not a full FEM restart; no interior eigenvectors or assembled matrices are persisted. The CLI is not yet wired to this API.


## 1.2.176 standalone synthesis SDK boundary

`bell_synthesis` is now an independently installable C++17 target exported as `BellSynthesis::bell_synthesis`. Its source list is restricted to profile/cache serialization, boundary-mode reconstruction, clapper contact, acoustic radiation/listener geometry and audio primitives. It does not link `bell_core`, FEM, optimizer, mesh generation, CalculiX or CLI sources. The external consumer example resolves the installed package and links only this target. This is a library boundary, not yet a complete real-time renderer: one-call rendering facade, streaming sink, stable C ABI and game-engine bindings remain planned.
