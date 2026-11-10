# Development plan

## Current baseline — 1.2.179

The modal continuation audit and geometric stereo comparison are complete within their documented scope. The source tree is consolidated; generated listening files and raw experiment matrices are not shipped. WAV serialization is centralized and tested. Mono and stereo far-field synthesis use the Fourier/Bessel directional-radiation representation. Version 1.2.170 added an adaptive local polar field family to the stereo path: only the local polar interval covering the supported head-yaw range (-90 to +90 degrees) is built, with maximum 2-degree node spacing, complex interpolation, and phase-only azimuth/frame rotation. FEM, clapper contact and modal time evolution remain shared between channels. The production optimizer and its objective remain unchanged.



## Next configuration milestone — generic acoustic target profiles

The next configuration change is a target-profile refactor, not a new cultural/genre-specific optimizer. Stage 8 still has a fixed four-frequency/four-amplitude target representation and hard-coded sorted-mode mappings. The target design is specified in `TARGET_PROFILE_ARCHITECTURE.md`: variable-length target terms, explicit reference selectors, metadata-only profile identity, and geometry seed selection owned by `evolution*.json`.

Implementation order:

1. Capture a deterministic baseline for the current production Western profile and preserve its score/ranking.
2. Introduce generic in-memory target vectors and a strict v2 JSON parser while retaining a v1 compatibility reader.
3. Move `seed_family` ownership to evolution configuration; treat legacy target `tradition` only as a controlled migration fallback.
4. Stage 10 implements generic variable-length v2 target vectors, v1 conversion, sorted-index selectors, explicit frequency/amplitude references, and decoupled `evolution.seed_family`. Remaining gates: baseline score/ranking parity, stronger parser path diagnostics, and degenerate-mode treatment.
5. Add variable-count, reference-selector, malformed-config, profile-name-independence and geometry/acoustic-separation tests.
6. Compare v1 and migrated v2 Western runs on identical candidates and settings before any default or production switch.

The first selector implementation may support only sorted mode indices. It must report that limitation and must not claim physical modal-family identification. Degenerate circular mode groups require a separate validated selector/evaluator path.

### Round-bell clapper integration — 1.2.179

Added `evaluateClapperModalStateFromCandidate(config, inner_profile_si, final_geometry_modes, result)`. This is the candidate-facing provider API: it accepts only precomputed modal modes and their matching SI inner-wall profile, so it cannot create a mesh or invoke FEM. The older `evaluateClapperModalStateForGeometry` remains a standalone diagnostic wrapper that owns its own FEM solve and must not be called from candidate evaluation. A dedicated contract test checks insufficient modes, invalid profile and output reset. The clapper objective now groups adjacent modes within a relative frequency gap of 5e-5, averages group frequency and combines pressure magnitudes by root-sum-square. Its contract test verifies invariance to a rotation of a degenerate two-mode basis. The full configured Release build reached 100%; the focused 11-test regression selection passed 11/11. The 97-test inventory was not run in full for 1.2.179.

These are still primitive/API changes only. The evolutionary optimizer is not yet wired to the clapper objective; `western_five_partial` remains the only production objective. The next step is to establish the physical mapping of circular FEM modes into the five target partial families, then make candidate evaluation obtain a final-physically-scaled circular FEM state, pass that exact state to the new API, and add an end-to-end candidate test.

### Round-bell clapper optimizer readiness — 1.2.177

The active optimizer scope is circular bells only. The standard configuration keeps `optimize_ellipticity=false`, pins `ellipticity=1.0`, and leaves `BUILD_CHINESE_RESEARCH=OFF`. The detached Chinese-bell research tree is frozen; no Chinese/elliptical objective is part of this release.

A readiness audit and regression pass were performed without changing the production objective. The complete configured Release build reached 100%. A 22-test round/clapper-focused selection passed 22/22, including modal tracking, shape/scale, mesh convergence, contact geometry/dynamics, modal evaluator, clapper objective primitive and end-to-end clapper WAV smoke. A further 12 integration/library contracts passed 12/12. The 96-test CTest inventory includes two long elliptical-convergence audits which are outside this round-only release gate; a broad non-elliptical suite run was interrupted by the execution time limit, so no claim of a complete 96/96 run is made.

The audit found two hard integration requirements before switching objective paths:

1. `evaluateClapperModalStateForGeometry` currently constructs a new cylindrical FEM model and solves it internally. Production candidate evaluation must instead pass the candidate's already-computed final-geometry `ModalData` and matching inner-wall profile to a reusable clapper evaluator, avoiding a duplicate FEM solve.
2. The current optimizer applies analytic scale-after-solve to modal frequencies and geometry. This is sufficient for the legacy dimensionless modal-ratio objective, but not for dimensional clapper contact. The clapper objective must receive modal fields recomputed on the final scaled geometry (or a separately proven dimensional transformation); simply reusing pre-scale fields with unchanged clapper dimensions is invalid.

The existing `western_five_partial` production score is unchanged. `evaluateClapperAcousticObjective` remains a tested primitive, not a production candidate objective. The demo clapper values remain uncalibrated. Before objective promotion, add a candidate-evaluation integration contract, final-scale FEM reuse contract, mesh/mode-count convergence, rank-stability comparison on a fixed candidate set, and explicit reporting of clapper parameter provenance.

The 1.2.171 real-FEM contract covers three geometry/mesh cases, 54 FEM modes spanning 79.55–735.10 Hz, 1,620 complex ear/direction comparisons, and 162 direct-quadrature spot checks. Worst per-mode peak-normalized complex error is 0.01182%; maximum phase error above 10% of each mode peak is 0.00127 degrees. Both pass the declared 1% / 1-degree gate by a wide margin. This is numerical validation of the representation/interpolator on the tested cases, not physical validation of the far-field model.

## Next work, in order

### 1. Stereo phase and timing — numerical contract complete

- The shared modal pressure kernel applies the centre-reference-to-ear propagation delay exactly once as a time shift.
- Geometric interaural time difference, onset causality and sample-grid quantization are covered by `stereo_phase_delay_contract`.
- Complex modal radiation phase is retained in the transfer coefficient; no second propagation phase multiplier is used.
- Keep far-field as default and finite-distance as an opt-in approximate comparison path until independent physical validation exists.

### 2. Perceptual/physical validation

- Compare with measured recordings only when bell geometry, material, clapper mass properties, strike conditions and microphone positions are known.
- Keep listening preference separate from physical accuracy and calibrated SPL.

### 3. Round-bell clapper optimizer integration — active next stage

The next optimizer work is restricted to **circular/round profiles**. The production geometry must use `ellipticity = 1` exactly; `optimize_ellipticity` remains false and `use_direct_elliptical_fem` remains false. The clapper-driven objective is to be validated and integrated as a separate selectable objective for round bells only. It must not silently replace the existing `western_five_partial` objective or change its scores.

Required sequence:

1. Establish a reproducible circular baseline with a fixed seed, material preset, target profile and solver settings.
2. Reuse the candidate's final-geometry FEM modes and inner-wall boundary fields for clapper excitation; avoid an independent duplicate FEM solve. The clapper calculation must run on the **physically scaled final geometry**, because its mass, inertia, tip radius and contact are dimensional quantities.
3. Validate contact integration and modal release on round profiles: time-step convergence, contact-event convergence, modal-frequency convergence, clapper-projected modal-amplitude convergence and sensitivity to retained mode count.
4. Verify modal-family identity and reject/flag ambiguous assignments rather than treating frequency rank as a guaranteed physical family. The 1.2.179 degenerate-group aggregation protects against arbitrary basis rotations inside a pair but does not solve family crossings or target-family identification.
5. Compare old and clapper objectives on the same fixed candidate set. Report both scores and rankings; do not change the default objective until the new path passes the declared gates and the user approves promotion.
6. Record the clapper parameter set as illustrative/uncalibrated unless it is supplied from measurement or CAD. A demo configuration is suitable for software contracts, not a claim of physical calibration.

Initial numerical acceptance gates for the round-bell pilot (to be reported per fixture, not averaged away):

- Halving the contact integration step changes peak contact force by less than 1% and contact time by less than 0.5%.
- Refining the bell mesh changes each target frequency ratio by less than 1% and each target clapper-driven modal-pressure ratio by less than 2%.
- Increasing retained modal count changes the objective score by less than 2%; otherwise the score is marked mode-count limited.
- Candidate ranking between the baseline and refined settings has Spearman correlation at least 0.95; every target mode must have an unambiguous tracked family assignment.
- The legacy `western_five_partial` score on the same candidate/solver inputs must remain numerically unchanged when the clapper objective is enabled as an alternative.
- The clapper simulation is run only after final physical scaling and must record the exact geometry, material preset, solver settings, clapper parameters, strike fraction and deterministic seed.

These are numerical/software gates, not proof of agreement with a real bell. Physical calibration still requires measured clapper and bell data.

The circular path must be independently testable with no Chinese-research target, elliptical geometry, stereo listener pose or repeated-strike synthesis involved.

## Chinese-bell research reopening condition

The current reduced model remains `RESEARCH_BLOCKED`. Reopen only with independent modal-shape/geometry evidence or a clearly falsifiable model change. Do not run another unconstrained sweep or tune modal tracking/contact width to manufacture the desired interval.

## Explicit non-goals

No 3-D FEM, HRTF/pinna, room acoustics, arbitrary alloy-property editor, automatic production-objective switch, Chinese-bell optimization, elliptical-profile development or repeated-strike/percussion synthesis in the round-bell optimizer stage. Chinese-bell work remains isolated under `research/` and is not part of the active development path.


### Stereo and reusable-profile gates completed through 1.2.174; SDK packaging added in 1.2.176

1. Real-FEM polar interpolation validation across profile/mesh cases: complete; measured errors remain below the declared 1% complex / 1-degree phase limits.
2. Reusable prepared receiver-query API with yaw/azimuth changes and stable polar-field build count: complete.
3. Geometry/material/solver-bound full-polar acoustic field cache and receiver preparation from cache: initial library API complete; real-FEM interpolation quality is covered by the modal-field contract. CLI cache flags and automatic cache management remain deferred until real workload size/build-time measurements.
4. Propagation-delay accounting, ITD, onset causality and once-only time-shift contract: complete.
5. Next: measure cache file size/build time on larger real FEM modal sets, then decide whether CLI cache-in/out is worthwhile. After that, revisit perceptual listening comparisons and physical measurement planning.


### 4. Acoustic/modal profile cache — library bundle API complete

- `AcousticFieldCache` stores geometry/material/solver-dependent modal radiation coefficients, not strike/clapper excitation or release state.
- `BellModalProfileBundle` pairs `.bmodal` structural boundary traces with `.acfield` directional radiation fields. The structural half includes inner/outer profiles, harmonic labels, Fourier radial/azimuthal/axial boundary traces, modal frequencies and modal masses.
- Real-FEM round-trip validation covers four actual modes from the 17-node reference geometry, preserving inner/outer boundary traces and loading a receiver-query state without rebuilding radiation surface integrals. The broader polar interpolation gate still covers 54 modes across three geometry/mesh cases.
- Cache identity must bind canonical geometry, material preset and FEM/modal solver settings.
- Full polar grid permits reuse across listener distance/direction/yaw; no surface integrations are repeated when preparing a receiver from a loaded cache.
- The bundle is a boundary-modal synthesis profile, not a full FEM restart: it does not serialize stiffness/mass matrices or interior eigenvectors. It is sufficient for the present clapper boundary-contact and radiation paths when the caller supplies the matching geometry/material/solver identity.
- The two files are a logical pair, not one atomic filesystem transaction; a missing half invalidates the profile.
- Next: measure file size/build time on larger real FEM modal sets, add cache-quality reporting and integrate explicit profile load/save into a fast synthesis entry point. Do not enable an implicit cache or guess fingerprints.


## Standalone synthesis SDK (1.2.176)

The first independent library boundary is in place: `BellSynthesis::bell_synthesis` can be installed and consumed without optimizer/FEM/CLI linkage. Next gates: (1) add a high-level profile-backed strike/render session, (2) add caller-owned sample-buffer rendering suitable for game audio callbacks, (3) add a C ABI only after the C++ API stabilizes, and (4) add a minimal external-host integration test that renders a short buffer from a fixture profile. Keep optimization and profile generation offline.


### Immediate round-bell clapper-objective blocker (1.2.180)

The combined-harmonic direct FEM solve does not currently provide a reliable one-label-per-mode family identity: the spectrum audit found substantial mixed Fourier mass in all 16 close-frequency candidate pairs. The next implementation must isolate azimuthal harmonic blocks while preserving axis regularity (`m=0`, `m=1`, and `m>=2`), attach explicit harmonic metadata, and track radial/axial families within each block. Only after that can the five target partials be mapped safely. The existing frequency-only degeneracy grouping must not be interpreted as proof of physical family identity. The production objective remains unchanged.


## Multi-fidelity optimizer investigation (round bells)

A first calibration probe is available in `MULTIFIDELITY_CALIBRATION.md`. It
compares the fast axisymmetric solver and direct cylindrical vector FEM on a
five-geometry sample. The low observed rank correlation (rho=0.30) and the
large frequency-rank mismatch show that the two solvers cannot yet be used as
coarse/fine evaluators with a simple sorted-frequency correspondence. The next
step is explicit modal-family correspondence and a larger deterministic
calibration set; only then may a conservative shortlist coordinator be added.
The production optimizer and `western_five_partial` objective remain unchanged.


## Multi-fidelity optimizer — current calibration state (1.2.182)

The local 30-geometry probe now associates fast and direct-vector-FEM `m=0`
modes using outer-wall displacement traces. This is a diagnostic association,
not yet a validated multi-fidelity ranking model: near-degenerate families
need subspace comparison, the geometry sample must cover the intended search
domain, and shortlist recall must be measured. Keep the production optimizer
on `western_five_partial`; do not enable fast-stage candidate rejection until
these gates have been met.


The 1.2.182 30-case probe flagged all 30 cross-solver assignments as ambiguous
under the provisional minimum trace-similarity threshold (0.75). The all-case
provisional Spearman value was 0.679 and the worst matched-frequency error was
102.66%; neither is a promotion result. The fast evaluator must not reject
candidates based on these associations. Next implement a degenerate-subspace
comparison for boundary traces and validate it against full-volume overlap;
otherwise retain conservative fallback to accurate evaluation.


## 1.2.183 — modal subspace validation primitive

Implemented a small, solver-independent rank-1/rank-2 principal-angle metric for
outer boundary traces. Synthetic tests verify basis-rotation invariance, rejection
of an orthogonal subspace, and rejection of rank-deficient bases. This does not
resolve the 1.2.182 cross-solver calibration failure: the metric is not yet wired
into assignment and has not yet been validated against full-volume mass-weighted
overlap. Keep the multi-fidelity shortlist disabled and retain exact FEM evaluation
for all production candidates until physical group identification, cross-metric
validation, and shortlist-recall gates pass.


## 1.2.184 — principal-angle invariance regression

The rank-2 principal-angle contract now includes 30 deterministic basis-change
cases (rotation, permutation, sign and scale) and 30 analytically calibrated
principal-angle perturbations. Both sweeps passed with worst absolute error
2.22e-16 in the GCC 14.2 Release run. The existing controlled cross-mesh FEM
fixture also passed with minimum one-to-one overlap 0.999854 and pair-subspace
similarity 0.999966. These results validate the numerical primitive and its
basis invariance; they do not validate cross-solver family assignment or the
fast/accurate ranking model. The 1.2.182 calibration result remains a blocker.
Keep candidate rejection/shortlisting disabled and run full accurate evaluation
for production candidates. Next gate: apply subspace evidence in the diagnostic
cross-solver assignment, compare it with full-volume mass-weighted subspaces,
and measure assignment stability and shortlist recall on a representative
round-bell geometry set before considering any production integration.


## Stage 6 update — round-bell cross-solver association (2026-11-10)

A research-only full-volume meridional matching probe now assigns the first five legacy modes to distinct candidates among the first eight direct cylindrical `m=0` modes. The same-mesh 61-geometry run produced 305 matches with minimum diagonal-lumped-mass overlap 0.999999999999 and maximum relative frequency difference 8.29e-13. Ten-geometry checks at 23 and 31 nodes also passed with overlaps at or above 0.999999999999. This points to a limitation in the previous all-eight boundary-trace assignment rather than a demonstrated mismatch in the first five meridional modes.

Do not promote this result yet: validate with the consistent FEM mass bilinear form, compare rank-2 subspaces for near-degenerate pairs, perform the boundary-trace square-vs-rectangular ablation, and expand refined-mesh coverage. `western_five_partial` remains the sole production objective; no optimizer candidate screening is enabled.
