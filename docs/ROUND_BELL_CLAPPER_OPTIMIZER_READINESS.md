# Round-bell clapper optimizer readiness — 1.2.179

## Scope

This audit is limited to circular/round bell profiles. The production optimizer configuration must keep `optimize_ellipticity=false` and `ellipticity=1.0`. `BUILD_CHINESE_RESEARCH` remains `OFF`. The Chinese-bell/elliptical research tree is detached and frozen; it is not a target of this work.

No production objective was switched. The production path remains `western_five_partial`. The clapper acoustic objective is a separately tested primitive and is not called by production candidate evaluation.

## Verification performed

- CMake Release configuration: GCC 14.2, C++17, `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF` — passed.
- Full configured Release build in `/mnt/data/bell_178_work/build_179` — all targets reached 100%.
- Focused CTest selection — 11/11 passed: modal-family tracking, optimizer shape/scale, cylindrical degenerate-subspace tracking, subspace similarity/conditioning, clapper dynamics, modal evaluator, candidate-data reuse, clapper acoustic objective and `internal_clapper_audio_smoke`.
- The 97-test inventory was not run in full for 1.2.179. The previous 1.2.178 broad run completed 95 tests successfully; two ellipse-specific convergence audits exceeded the execution window. No new full-suite pass is claimed.

## Code-path findings

### Existing production objective

`EvolutionaryOptimizer::computeError` still scores the legacy five-partial target using frequency ratios, excitation-amplitude ratios and a geometry proportion penalty. The frequency target is reached by analytic scaling after the modal solve. The amplitude term has solver-dependent semantics and is not a clapper-driven response.

### Clapper primitive

`evaluateClapperAcousticObjective` consumes a `ClapperModalState` and scores modal-family pressure ratios. Version 1.2.179 groups adjacent frequencies within a relative tolerance of 5e-5, averages their frequencies, and combines their pressure magnitudes by root-sum-square. This makes an exactly degenerate circular-bell subspace invariant to orthogonal changes of eigenvector basis. The contract test now checks a two-mode degenerate Hum family, a rotated response basis, common-gain invariance, and rejection of a zero target-family response. This still does not establish candidate-level optimizer integration, calibrated acoustic accuracy or robust ranking.

### Duplicate FEM solve

`evaluateClapperModalStateForGeometry` currently calls `make_cylindrical_bell_model` and `evaluate_cylindrical_bell_modes` itself. If called from `evaluateCandidate` as-is, it would solve FEM a second time for the same candidate. The integration API should be refactored to accept:

- the candidate's final `ModalData`, including the modal boundary fields required by contact;
- the matching SI inner-wall profile;
- the typed clapper/contact configuration.

The geometry-level convenience wrapper can remain for standalone synthesis/tests, but production optimization should use the precomputed-data overload.

### Physical scaling boundary

The current optimizer's analytic scale-after-solve is suitable for the legacy dimensionless frequency-ratio objective. It is not automatically valid for the dimensional clapper system. Mass, inertia, pivot geometry, tip radius and contact compliance are dimensional. The clapper objective must use FEM modes and inner-wall fields computed on the final physically scaled geometry. Reusing pre-scale modes while keeping clapper dimensions unchanged is not acceptable.

### Modal identity and degenerate families

The 1.2.176 tracking-order correction is in place. Version 1.2.179 makes the clapper objective invariant to basis rotation within a numerically degenerate adjacent-frequency group. This is a local subspace aggregation rule, not a full candidate-to-candidate family tracker. It does not resolve accidental crossings between distinct families, nor does it prove that the first five frequency groups correspond to the intended five acoustic partials for the production FEM basis. Candidate-level scoring must still reject or explicitly flag ambiguous family assignment; frequency rank alone must not be treated as a stable physical family label.

## Physical calibration boundary

The bundled 5 kg / 0.46 m COM distance / 0.52 m tip-centre distance / 35 mm tip radius / 3.2 rad/s initial angular velocity set is demonstrative and uncalibrated. It is suitable for deterministic software tests and sensitivity studies only. A physical prediction needs actual clapper mass properties, pivot and tip geometry, release/impact state, bell material preset and measured geometry.

The current spherical Hertz law assumes a locally flat wall and elastic, undamped contact. It omits wall-curvature correction, friction, plasticity and contact damping. These are model limitations, not defects that can be fixed by tuning the optimizer weights.

## Required work before enabling the new objective

1. **Completed in 1.2.178:** add a provider entry point consuming already-computed modes and the matching inner profile, with no geometry/material/FEM dependency; add an API contract test. The standalone geometry-level wrapper remains for diagnostic synthesis.
2. **Completed in 1.2.179:** make the clapper acoustic objective invariant to basis rotations inside numerically degenerate adjacent-frequency groups using average group frequency and root-sum-square pressure. This does not itself establish physical family labels.
3. Establish and test the mapping from the circular FEM spectrum to the intended five acoustic partial families, including degenerate pairs and near crossings.
4. Wire candidate evaluation to compute or obtain modal fields on the final physically scaled circular geometry, then pass those exact modes/profile to `evaluateClapperModalStateFromCandidate`. Do not use analytically frequency-scaled pre-scale modes for dimensional contact.
5. Add an end-to-end candidate-evaluation test with deterministic seed and explicit `objective_model=clapper_acoustic` selection; preserve `western_five_partial` as the default.
6. Validate contact time-step convergence: peak force change <1%, contact time change <0.5% when halving the step.
7. Validate mesh convergence: target frequency ratios <1% change and clapper-driven modal-pressure ratios <2% change.
8. Validate retained-mode convergence: objective score <2% change, otherwise mark the candidate mode-count limited.
9. Compare candidate ranking at baseline/refined settings; require Spearman rho >=0.95 and no ambiguous target-family assignment.
10. Record objective name, physical geometry after scaling, solver/mesh/mode count, material preset, clapper parameters and provenance, strike fraction, time step and deterministic seed in every run report.
11. Compare old and new scores/rankings on the same fixed candidate set. Do not promote or replace the production objective until numerical gates pass and physical target data is defined.

## Release conclusion

The round-bell FEM, modal, contact and acoustic building blocks are in a sufficiently healthy state to continue integration. Version 1.2.178 closes the provider API's duplicate-FEM hazard, and 1.2.179 makes degenerate-family pressure scoring basis-invariant. The 1.2.179 configured Release build reached 100%, and the focused regression selection passed 11/11. The previous broad 1.2.178 run completed 95/97 tests; two ellipse-specific convergence audits exceeded the execution window. The optimizer does not yet call the new entry point. Degenerate modal subspace scoring is now basis-invariant at the objective primitive, but target-family mapping, final-scale modal recomputation, candidate-level scoring, ranking validation and physical calibration remain open. The new objective is **not ready to become the production default**.


## 1.2.180 modal-family spectrum audit

The clean 1.2.179 build exposed a gap in the previous evidence trail: the earlier `build_179` directory was configured against 1.2.178. A fresh build was therefore configured directly from 1.2.179 sources and the focused clapper/modal tests were rerun. The new 1.2.180 diagnostic test uses the actual circular reference geometry and direct cylindrical FEM. It finds 16 close-frequency candidate pairs, all flagged as mixed/ambiguous by a conservative dominant-harmonic energy threshold. This means neither frequency rank nor the largest Fourier component is currently a defensible physical family label. The frequency-only degeneracy aggregation from 1.2.179 is a local basis-invariance improvement, not a complete family identity solution.

Next: generalize assembly and axis regularity so each azimuthal harmonic block can be solved independently; preserve explicit harmonic-order labels, validate the `m=0`, `m=1`, and `m>=2` axis constraints, then perform within-block family tracking. Do not connect the clapper objective to production candidate ranking before this and the final-scale FEM/calibration gates are complete.


### 1.2.180 build verification

A fresh Release build was configured from the 1.2.180 source directory with GCC 14.2, C++17, `BUILD_SMOKE_TESTS=ON`, and `BUILD_CHINESE_RESEARCH=OFF`. All configured targets built successfully to 100%. The focused optimizer/modal/clapper selection passed 13/13, including the new real-geometry spectrum-family audit and the internal-clapper audio smoke test. This is focused evidence only; it is not a full-suite claim and does not close the physical calibration or candidate-ranking gates.
