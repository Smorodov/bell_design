# Changelog

## 1.2.184 — modal subspace invariance audit

- Expanded the principal-angle contract with two deterministic 30-case sweeps: arbitrary in-plane basis rotations/permutations/sign changes/scalings, and controlled out-of-plane principal-angle perturbations.
- The 30 basis-invariance cases retained unity similarity with worst absolute error 2.22e-16; the 30 calibrated-angle cases matched the analytical cosine with worst absolute error 2.22e-16.
- Built and passed `modal_subspace_similarity_contract`, `modal_subspace_conditioning_contract`, `modal_tracker_fem_continuation`, and `cylindrical_modal_field_cross_mesh` with GCC 14.2 / C++17 Release.
- Cross-mesh controlled FEM fixture retained minimum one-to-one overlap 0.999854 and two-mode subspace similarity 0.999966.
- This is a numerical contract expansion, not a resolution of the 1.2.182 cross-solver calibration failure: the trace metric is not yet used for cross-solver assignment, and no 30-geometry recall promotion is claimed. Multi-fidelity shortlist remains disabled; production objective remains `western_five_partial`.

## 1.2.183 — degenerate boundary-trace subspace metric

- Added a rank-1/rank-2 principal-angle similarity metric for boundary-trace subspaces.
- The metric is invariant to eigenvector sign, scaling, and orthogonal basis rotation; invalid/rank-deficient bases are rejected.
- Added synthetic tests for rotation invariance, orthogonal subspaces, and rank deficiency.
- This is a diagnostic primitive only: cross-solver family assignment and optimizer screening remain disabled pending validation against full-volume mass-weighted overlap and shortlist recall.
- Production objective remains `western_five_partial`; Chinese research remains excluded from the normal build.

## 1.2.182 — cross-solver modal-shape calibration probe

- Expanded the local round-bell multi-fidelity probe from five to 30 deterministic profile perturbations.
- Enabled outer-wall modal-shape capture in both the fast meridional solver and direct cylindrical `m=0` vector FEM reference.
- Associated modes by sign-invariant outer-wall radial/axial trace similarity using the existing `trackModalFamilies` assignment; frequency is only a weak tie-breaker.
- Reports matched-family frequency errors, trace similarity, ambiguous assignments and Spearman ranking correlation.
- Explicitly documents that boundary-trace matching is not a full-volume or degenerate-subspace proof; no production shortlist or objective switch is enabled.


## 1.2.180 — round-bell modal-family spectrum audit

- Added `clapper_modal_spectrum_family_audit`, a reproducible direct vector-FEM diagnostic for the circular reference bell.
- The audit shows that a dominant Fourier label is not a reliable family identity in the current combined-harmonic eigensolve: all 16 close-frequency candidate pairs in the 40-mode diagnostic were flagged as mixed/ambiguous by the 0.95 mass-energy criterion.
- Documented why frequency-only degeneracy grouping is not sufficient to prove physical modal-family identity or the five-partial mapping.
- Defined the next required step: independently solve harmonic blocks while preserving the correct axis constraints, then track families within each block.
- No production optimizer objective was changed; `western_five_partial` remains the default. The clapper objective remains experimental.

# 1.2.179 — degenerate modal-family objective

- The clapper acoustic objective now groups adjacent modes whose relative frequency gap is at most 5e-5.
- Group frequency is the arithmetic mean; group pressure is the root-sum-square of member modal-pressure magnitudes, making an exactly degenerate subspace score invariant to orthogonal eigenvector-basis rotations.
- Updated the objective contract test to cover a degenerate Hum pair, basis rotation, common gain and zero target-family rejection.
- Updated round-bell readiness documentation: the grouping primitive is not a substitute for explicit physical target-family mapping or candidate-level modal tracking.
- Production optimizer remains on `western_five_partial`; no candidate-evaluation integration or default-objective switch in this release.


## 1.2.178 — reuse candidate modal state for clapper evaluation

- Added `evaluateClapperModalStateFromCandidate`, which consumes the already-computed modal modes and matching SI inner-wall profile. It has no mesh/model/material/FEM arguments and cannot launch a duplicate FEM solve.
- Retained `evaluateClapperModalStateForGeometry` as a standalone diagnostic convenience wrapper; documented that it must not be called from optimizer candidate evaluation.
- Added `clapper_candidate_data_reuse_contract` for insufficient mode count, invalid profile and output-reset behavior.
- Updated round-bell readiness and development-plan documents.
- Full configured Release build reached 100%. Of 97 CTest cases, 95 passed; two long ellipse-specific convergence audits exceeded the execution window. All five clapper-focused tests, including the new API contract and internal-clapper WAV smoke, passed.
- The production optimizer remains on `western_five_partial`; candidate wiring, final-scale FEM recomputation, ranking validation and physical calibration are still pending.


## 1.2.177 — round-bell optimizer readiness audit

- Froze active development scope to circular/round bell profiles; the detached Chinese-bell research module remains disabled by default and is not advanced by this release.
- Added a focused readiness report covering round-bell FEM/modal state, clapper contact, objective boundary, physical scaling and required integration contracts.
- Verified a complete configured Release build with GCC 14.2 / C++17.
- Passed the 22-test round/clapper-focused regression selection and a further 12 integration/library contracts. The complete 96-test suite was not claimed: two long elliptical-convergence audits are outside the round-only gate, and a broad run exceeded the execution window.
- Confirmed that production optimization still uses `western_five_partial`; the clapper acoustic objective remains unconnected.
- Identified required next implementation changes: reuse candidate FEM modal data rather than solving twice, and run clapper contact against modal fields from final physically scaled geometry.
- No objective formula, FEM equation, contact law or default optimization behavior was changed.


## 1.2.176 — optimizer modal-order audit and correction

- Audited the production objective boundary and documented the exact frequency, amplitude and geometry terms, including the solver-dependent physical meaning of `ModalMode::amplitude`.
- Fixed modal-family tracking order: tracking now runs after direct elliptical FEM evaluation or coordinate morphing, so later replacement/sorting no longer silently discards the tracked family assignment before objective scoring.
- Added `docs/OPTIMIZER_AUDIT.md` with remaining physical, convergence, acoustic-reference and acceptance gates for future clapper-driven optimization.
- Production objective remains `western_five_partial`; the clapper objective remains a tested but unconnected primitive. No objective switch or Chinese objective integration was made.


## 1.2.172 — prepared stereo receiver state

- Added immutable `PreparedStereoRadiation` to separate per-mode polar-field construction from listener-yaw queries. It returns complex left/right transfer, ear range scales and propagation delays, and rejects yaw outside its prepared envelope.
- Routed production clapper far-field stereo synthesis through the prepared receiver object after FEM/contact/modal release. Mono and finite-distance paths remain unchanged.
- Extended the real-FEM polar contract with yaw/azimuth query matrix across three geometry/mesh cases, exact-polar comparisons, stable polar-field build-count assertions and out-of-envelope rejection.
- Full Release build passed; focused stereo/radiation tests passed 9/9; 92 CTest cases registered, full suite not run. Production stereo WAV smoke passed (48 kHz, 16-bit, two channels, 32 modes).
- No FEM/contact physical equations or optimizer objective changed. No interactive UI or persistent cross-process cache is claimed.


- Added an adaptive local polar family of complex directional radiation fields covering the supported head-yaw interval (-90 to +90 degrees), with maximum 2-degree node spacing and complex interpolation.
- Stereo far-field synthesis now evaluates each ear at its own polar and azimuth angles without repeating FEM, clapper contact, or modal time integration. The local grid typically needs 2–4 radiation surface integrations per mode and can be reused across supported head yaw. Mono remains on the single-polar path; finite-distance mode is unchanged.
- Added interpolation and frame-rotation contract; synthetic fixture error is 0.442% of reference peak on the global 91-node grid and 0.365% across the production yaw envelope.

## 1.2.169 — clapper/modal covariance and polar-angle audit

- Added an end-to-end clapper frame-covariance contract using the production contact integrator, an exactly degenerate cosine/sine modal pair, a rotated strike azimuth, release modal state and complex directional radiation.
- Added a stereo polar-angle approximation audit. In its documented synthetic fixture, the centre-polar shortcut error is 29.30% at 1.5 m and 14.20% at 3 m; a polar-specific field agrees with direct quadrature to numerical precision.
- Updated the acoustic audit, validation record and development plan. The current production shortcut remains unchanged and is explicitly treated as an approximation pending a 3D/polar-dependent field.
- No production FEM, clapper, acoustic-kernel or optimizer-objective equations changed.

# Changelog

## 1.2.183 — degenerate boundary-trace subspace metric

- Added a rank-1/rank-2 principal-angle similarity metric for boundary-trace subspaces.
- The metric is invariant to eigenvector sign, scaling, and orthogonal basis rotation; invalid/rank-deficient bases are rejected.
- Added synthetic tests for rotation invariance, orthogonal subspaces, and rank deficiency.
- This is a diagnostic primitive only: cross-solver family assignment and optimizer screening remain disabled pending validation against full-volume mass-weighted overlap and shortlist recall.
- Production objective remains `western_five_partial`; Chinese research remains excluded from the normal build.

## 1.2.166 — shared audio output and far-field radiation

- Extracted mono/stereo PCM16 RIFF/WAVE writing into `src/audio/wav_writer.*`; all synthesis workflows now share it. Added validation for headers, interleaving, clipping, invalid inputs and file-write failures.
- Moved direct azimuthal radiation and Fourier/Bessel field construction into `src/acoustics/directional_radiation.*`. Production mono and stereo far-field synthesis now share the harmonic field builder.
- Added a complex-amplitude/phase equivalence test. Across the tested frequency/polar/azimuth matrix, maximum relative complex difference was `2.99e-14`. The direct integral remains as a regression reference only.
- Full configured build succeeded; focused 10/10 WAV/acoustic/contact/stereo tests passed. Full suite and Windows/MSVC remain untested.
- No FEM/contact physical formulas or optimizer objective changed.

## 1.2.165 — first source-level cleanup pass

- Enabled the same strict warning level for the CLI executable as for the core libraries (`-Wall -Wextra -Wpedantic`, or `/W4` on MSVC).
- Removed two proven-unused internal helpers: `basisForLabel` in contact synthesis and `smoothstep` in the legacy synthesis module.
- No physical formulas, acoustic transfer kernels, clapper/contact behavior or optimizer objective changed.
- Rebuilt all configured targets without compiler warnings and reran the focused 10-test regression group.


This is a concise milestone history, not an append-only list of every patch release.

## 1.2.164 — source and evidence consolidation

- Removed generated WAVs, raw CSV matrices, metadata sidecars and historical per-version reports from the source package.
- Replaced research outputs with maintained numerical summaries and explicit interpretation/limitations.
- Consolidated architecture, modal tracking, clapper/contact, stereo and validation documentation.
- Removed three one-shot result-generating audit tests after preserving their findings.
- No production FEM, contact, optimizer objective or acoustic kernel formula changed.

## 1.2.161–1.2.163 — stereo comparison studies

- Added real-FEM far-field/finite-distance receiver comparison, matched stereo WAV comparison and a 27-configuration listener geometry matrix.
- Established shared FEM/contact/modal history for left and right channels.
- Retained far-field as default; finite-distance remains opt-in because neither approximate kernel has independent physical validation.

## 1.2.159–1.2.160 — modal continuation and acoustic mesh study

- Added controlled one-to-one cross-mesh modal assignment and principal-subspace comparison.
- Closed the diagnostic modal-tracking stage for the tested geometry, q path, mesh and Fourier cases.
- Quantified selected modal-pair acoustic mesh refinement; coarse mesh is not accepted as converged for that metric.

## 1.2.140–1.2.158 — clapper and acoustic foundation

- Added inner-wall modal export, spherical-tip contact geometry, rigid-pendulum clapper dynamics, conservative Hertz contact and post-separation modal audio synthesis.
- Added clapper-to-modal evaluation and a clapper acoustic objective primitive, without switching the production optimizer.
- Added frequency-dependent modal decay reference, directional acoustic field, stereo receiver geometry and opt-in finite-distance radiation comparison.

## Earlier solver foundation

Earlier releases introduced the internal sparse eigensolver, cylindrical Q4 vector FEM, Fourier modal families, exact cylindrical-axis regularity, localized strike projection and mesh/modal validation. Detailed equations and contracts are in `docs/MATHEMATICAL_AND_ALGORITHMIC_FOUNDATIONS.md`, `docs/MODAL_ANALYSIS.md` and `docs/VALIDATION.md`.


## 1.2.167 — coordinate-covariant directional-field reuse
- Added no-copy rotated-frame evaluation for immutable-base directional acoustic fields.
- Expanded rotation contracts with joint-frame covariance, additive composition, repeated receiver-angle equivalence, and base-field immutability checks.
- Documented the boundary between algebraic covariance and end-to-end production geometry validation.

## 1.2.168 — joint stereo coordinate covariance
- Added a contract for common source/strike/listener z-axis rotation, complex left/right transfer invariance, ear geometry and propagation-delay invariance.
- Documented the remaining production clapper/modal-phase integration gap and far-field ear-polar-angle approximation.
- Full configured build succeeded; 86 of 88 registered tests are confirmed passing. Two long ellipticity audits exceeded the available test-run window and remain unverified.


## 1.2.171 — real-FEM polar interpolation validation

- Added `real_fem_polar_interpolation_contract`, validating local 2-degree complex polar interpolation using actual cylindrical FEM modal fields across a 17-node baseline, a 25-node refinement, and a second bell profile.
- Covered 54 modes, 1,620 complex ear/yaw/distance comparisons, and 162 direct-quadrature spot checks. Worst per-mode peak-normalized complex error was 0.01182%; maximum phase error above 10% of each mode peak was 0.00127 degrees. Both pass the 1% / 1-degree acceptance gate.
- Retained linear complex interpolation and the 2-degree maximum node spacing. No optimizer objective, FEM, clapper contact or physical radiation equations were changed.
- Updated stereo audit, validation, development plan, codebase audit and README.

- Full configured Release build passed with GCC 14.2 / C++17. Nine focused stereo/radiation tests passed 9/9; 92 CTest tests are registered (the full suite was not run). A production 0.25 s stereo clapper smoke generated a valid 48 kHz, 16-bit, 2-channel PCM WAV using 32 modes.

## 1.2.173 — acoustic field cache and stereo phase/time kernel

- Added versioned, checksummed binary `AcousticFieldCache` storing full-polar complex modal radiation fields and required geometry/material/solver fingerprints.
- Cache intentionally excludes strike force, clapper/contact history, modal release amplitudes/velocities, decay/T60, receiver geometry and audio samples.
- Added receiver preparation from a loaded cache; listener distance, direction, ear spacing and yaw can vary without new radiation surface integrations.
- Factored modal pressure evaluation into a shared mono/stereo kernel; propagation delay is applied once as a time shift, preserving complex directional modal phase.
- Added cache round-trip/identity/corruption/reuse contract and stereo phase/ITD/onset contract.

## 1.2.174
- Added a reusable boundary-modal profile containing inner/outer geometry and Fourier modal traces, plus modal frequency/mass metadata.
- Added a paired `.bmodal` + `.acfield` profile bundle with identity checks, checksums, dimensional validation and reconstruction of boundary-only `ModalMode` records.
- Strike force, clapper motion/contact, modal release state, T60 overrides and PCM remain outside the reusable profile.
- Added `bell_modal_profile_contract` covering round-trip preservation, profile identity rejection and receiver preparation from loaded radiation data. The real-FEM polar contract also round-trips four actual FEM modes and checks exact inner/outer Fourier boundary-trace preservation.
- Full Release build passed; 14 focused cache/FEM/stereo/clapper tests passed 14/14. 95 CTest tests are registered; the full suite was not run. A 0.25 s stereo clapper smoke generated a valid 48 kHz PCM16 stereo WAV using 32 modes.


## 1.2.175 — standalone Bell Synthesis SDK

- Added installable C++17 `bell_synthesis` library target and exported `BellSynthesis::bell_synthesis` CMake package.
- SDK target contains profile/cache I/O, boundary-mode reconstruction, clapper/contact integration, directional/listener radiation helpers and audio primitives. It does not link optimizer, FEM, mesh generation, CalculiX or CLI code.
- Added public header installation, package version metadata and an external CMake consumer example.
- Added a link contract that links only the SDK target.
- No optimizer objective or physical model coefficients changed.
- A high-level one-call renderer, streaming callback, stable C ABI and game-engine bindings remain future work.


## 1.2.181 — round-bell multi-fidelity calibration probe

- Added a reproducible five-geometry diagnostic comparing the current fast axisymmetric solver with the direct cylindrical vector FEM `m=0` path.
- The first probe reports Spearman ranking correlation 0.30 and up to 43.16% difference among the first five frequency ranks. The result shows that raw frequency-rank mapping is not sufficient for multi-fidelity screening.
- Added `docs/MULTIFIDELITY_CALIBRATION.md` with setup, observed metrics, limitations, and next calibration requirements.
- No production optimizer behavior, objective, solver default, or Chinese/elliptical scope changed.


## Documentation review — generic acoustic target profiles (2026-11-10)

- Added the selector-based target-profile architecture and explicit legacy-v1 migration contract.
- Corrected user/developer documentation to state that Stage 8 still has four fixed frequency/amplitude target slots and hard-coded sorted-mode mappings.
- Clarified that acoustic target profiles and geometry seed families are conceptually independent, while Stage 8 retains a legacy `tradition`-to-seed coupling.
- Clarified that `tuning_model` is currently allow-listed metadata and must not select a hard-coded objective in the generic architecture.
- Preserved the Sysoi target as a single-exemplar experimental reference, not a universal Orthodox tuning law.
- Added `docs/ROUND_BELL_TARGET_PROFILE_ARCHITECTURE_REVIEW_STAGE9_2026-11-10.md` with the source findings and next implementation gates.
- No production source code, objective, default profile, or project version changed in this documentation-only review.

## Stage 10 — generic target profiles (2026-11-10)

- Implemented selector-based v2 target profile loading with variable-length frequency/amplitude targets and explicit reference modes.
- Preserved legacy v1 files through conversion to the generic target representation.
- Unified production candidate scoring and removed profile-name-dependent objective/reporting branches.
- Moved geometry seed selection to `evolution.seed_family`, with legacy fallback and conflict detection.
- Added v2 example and loader contract test. Full Western score/ranking parity, full CTest and Windows verification remain open gates.
