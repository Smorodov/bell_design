# Production optimizer audit — 1.2.177 (round-bell readiness)

## Scope and safety boundary

This audit concerns the production Western `western_five_partial` optimization path. It does not authorize replacing that objective with the internal-clapper objective, nor does it promote the experimental Chinese-bell objective. The production objective remains unchanged in this release.

## Current production score (legacy fixed-slot target contract)

`EvolutionaryOptimizer::computeError` combines three terms:

1. Four log-relative frequency-ratio errors for mode indices `{0, 2, 3, 4}` relative to mode index `1` (Prime).
2. Four log-relative amplitude-ratio errors for indices `{1, 2, 3, 4}` relative to index `0` (Hum).
3. A soft penalty for mouth-diameter/height ratio outside the configured admissible interval.

The target frequency `prime_freq_target` is enforced using sorted mode index `1` by analytic geometric scaling after the modal solve. The scale factor is `f_prime_raw / f_prime_target`; all modal frequencies are divided by this factor and all geometry lengths are multiplied by it. Therefore the frequency-ratio terms are shape objectives, not independent absolute-frequency objectives.

## Physical meaning of the amplitude term

For the internal axisymmetric modal solver, `ModalMode::amplitude` is normally populated from the absolute value of

`part_rz * acoustic_volume_displacement * omega`,

then normalized by the maximum across returned modes. This is an acoustic monopole/meridional-excitation proxy. It is not the response to a specified internal clapper: it contains no clapper mass/inertia, pivot geometry, initial swing, strike height, local contact footprint, contact law integration or post-contact modal release state.

The CalculiX compatibility parser cannot reconstruct the surface acoustic coefficient from its legacy result file and falls back to a mechanical participation proxy. Consequently the same amplitude field can have different physical semantics depending on the selected modal solver. The CLI defaults to the internal solver, but this distinction must remain explicit in reporting and validation.

The implemented `evaluateClapperAcousticObjective` is a separate primitive. It consumes a previously computed `ClapperModalState` and scores modal pressure magnitudes/quadratures, but it is not called by production candidate evaluation. Its existence and unit contract do not establish calibrated clapper-driven optimization.

## Modal identity finding and correction

Before 1.2.176, `trackModalFamilies` was called immediately after the initial solver result. In ellipticity-enabled paths, that ordering was incorrect:

- the direct elliptical FEM branch replaced the already-tracked `raw_data` with a new frequency-sorted result;
- the coordinate-morph fallback explicitly sorted the morphed modes by frequency after tracking.

Both operations silently discarded the tracked family order before the objective indexed modes by `{0..4}`. Version 1.2.176 moves family tracking after the selected elliptical transformation/morphing path has finalized the modal set. The first transformed solve establishes the reference for that run; subsequent candidates are reordered before target scoring.

This fixes ordering in the existing algorithm. It does not solve all modal-identification limits: the current cross-candidate matcher needs compatible boundary sampling and uses individual shape similarity. Near-degenerate subspace rotations and changes in discretization remain separate validation topics.

## Remaining blockers before a clapper-driven optimizer can be accepted

1. **Physical inputs:** measured or explicitly bounded clapper mass, pivot-to-COM distance, pivot inertia, tip-centre distance/radius, initial angle/velocity and strike sector.
2. **Contact validity:** quantify the locally-flat Hertz approximation, contact time-step convergence, wall-curvature/contact-patch error and the effect of omitted damping/friction/plasticity.
3. **Modal convergence:** converge target frequency ratios, clapper-projected modal amplitudes and the final radiated pressure spectrum over mesh density, retained mode count and time step.
4. **Acoustic reference:** define receiver geometry, radiation approximation, bandwidth, time window, decay and whether the target is relative pressure, sound exposure or calibrated SPL.
5. **Objective design:** specify which modal families are targets, how untracked/missing modes are penalized, how phase and decay enter, and how common gain is treated. Do not assume modal pressure magnitude alone captures perceived bell quality.
6. **Robustness:** test candidate ranking under modest changes in material preset, mesh, clapper inputs and strike location; report uncertainty rather than hiding it in one scalar score.
7. **Acceptance gates:** set numerical tolerances before enabling the new objective in production.

## Next audit actions

- Add strict validation for `TargetData` and `EvolutionParams` at the optimizer boundary (finite positive ratios/frequency, finite nonnegative weights, feasible population/elite/generation and geometry-search ranges).
- Add a production candidate-evaluation contract that checks mode-family order survives both direct elliptical FEM and coordinate morphing paths.
- Report amplitude semantics and solver path in optimization output/result metadata.
- Establish a reproducible baseline run with fixed seed and archived input configuration before changing any objective equation.

## Explicit non-goals

No production objective switch, no Chinese objective integration, no 3-D FEM, and no repeated-clapper-strike/percussion synthesis in this audit stage.


## Scope correction for the next optimizer stage

The active optimization target is the **round/circular bell**. In the production run configuration, `optimize_ellipticity` must remain false, so `ellipticity` is pinned to exactly 1.0 and the direct elliptical FEM branch is not entered. The detached Chinese-bell research module and its elliptical objective are not part of this work and must remain disabled by default (`BUILD_CHINESE_RESEARCH=OFF`).

The next clapper integration must therefore be evaluated on round profiles only. Its FEM/modal data must correspond to the final, physically scaled geometry before the dimensional clapper contact simulation is run. The current analytic scale-after-solve shortcut is valid for dimensionless modal ratios, but it does not make a dimensional clapper configuration scale invariant. A clapper objective that consumes pre-scale modal fields with unchanged clapper dimensions is not acceptable.

The existing `western_five_partial` objective remains the reference baseline. A clapper-driven round-bell objective must first be added as an explicitly selected alternative, with side-by-side candidate scores and rankings. No default switch is implied by this audit or by the existence of the clapper synthesis command.


### Proposed numerical gates for the round-bell pilot

- Halve the contact integration step: peak contact force changes by <1%, contact time by <0.5%.
- Refine the bell mesh: each target frequency ratio changes by <1%, each target clapper-driven modal-pressure ratio by <2%.
- Increase retained modal count: objective score changes by <2%; otherwise report the result as mode-count limited.
- Compare candidate ranking at baseline/refined settings: Spearman rank correlation >=0.95, with no ambiguous target-family assignments.
- Verify the legacy `western_five_partial` score remains numerically unchanged when the clapper objective is selected as an alternative.

These are numerical acceptance gates only. They do not establish physical accuracy without measured bell/clapper reference data.


## Acoustic target configuration audit — Stage 8

### Historical Stage 8 finding (superseded by Stage 10)

The Stage 8 CLI loads target JSON into `TargetData`, whose four-element `std::array` members impose a fixed target count. The JSON loader copies into those arrays; it is not a generic variable-length profile loader. The production score still uses frequency indices `{0,2,3,4}` relative to reference index `1` and amplitude indices `{1,2,3,4}` relative to reference index `0`. The `tuning_model` allow-list and Orthodox-specific output branches are metadata/CLI/reporting coupling that must be removed from objective selection. The target loader currently also reads `tradition` from target JSON to select the geometry seed provider, although the desired architecture places `seed_family` in evolution/search configuration.

The next target-profile implementation is specified in `docs/TARGET_PROFILE_ARCHITECTURE.md`. It must first preserve the legacy Western score numerically, then support variable-length target vectors and explicit reference/mode selectors. The initial selector capability may remain `sorted_index`; it must not claim validated mode-family identity. The generic target evaluator must be profile-name agnostic.

## Stage 10 update — generic target evaluator

The Stage 8 fixed-slot objective described above is now the legacy v1 compatibility mapping, not the only production representation. Stage 10 converts v1 profiles to variable-length target vectors and accepts v2 selector-based profiles. Production scoring now evaluates each configured frequency and amplitude target through one generic evaluator; absolute-frequency scaling uses the configured reference selector. The v2 implementation currently supports `sorted_index` only. The separate experimental `clapper_acoustic_objective` helper still has a fixed legacy target interface and is not part of production candidate scoring. See `TARGET_PROFILE_IMPLEMENTATION_STAGE10_2026-11-10.md` for verified scope and outstanding parity gates.
