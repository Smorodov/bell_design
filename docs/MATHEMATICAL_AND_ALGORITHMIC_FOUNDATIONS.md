# Mathematical and algorithmic foundations

**Applies to BellOptimizer 1.2.113.**

Sparse generalized eigenpairs are accepted only after checking the post-Jacobi Ritz vectors in the normalized mass metric. The residual is a normwise backward error for `Kx = lambda Mx` and is evaluated for near-zero eigenvalues as well; no eigenvalue-floor exemption is used in the final convergence test. The shifted PCG operator has the form `K + delta M`, so the shift preserves generalized eigenvectors and is used only to improve the linear solve. This document describes the implemented model and explicitly separates production behavior from diagnostic/experimental paths.

## 1. Model pipeline

`BellGeometry → meridional mesh → FEM K/M → eigenpairs → modal tracking → excitation → acoustic transfer → synthesis`.

The non-axisymmetric structural model uses a meridional Q4 mesh and a finite vector-valued Fourier expansion in azimuth. It is not a full 3-D production mesh.

## 2. Geometry and units

The bell profile is represented by a meridional outer contour and a normal-offset inner contour. Geometry validation checks positive wall thickness, offset validity and element quality. The geometry/FEM adapter converts the application's millimetre-based geometry and material values to SI before assembly.

Ellipticity is `q=b/a`, with `0<q≤1` and `q=1` circular. The mapped radius is

`r(φ,z) = α(φ) r_ref(z)`,

`α(φ) = 1 / sqrt(cos²φ + sin²φ/q²)`.

The mapping affects the element geometry/Jacobian and therefore both stiffness and mass matrices. It is not a frequency-only correction.

## 3. Vector cylindrical FEM

The displacement field is `u=(u_r,u_φ,u_z)`, expanded in the selected Fourier harmonics with sine/cosine channels. Odd harmonics are supported; in particular, `n=1` terms are required to represent global translations in a cylindrical vector basis.

The element matrices follow

`K_e = ∫ Bᵀ D B dV`,

`M_e = ∫ ρ Nᵀ N dV`,

with the physical cylindrical volume Jacobian. The global sparse matrices are assembled from element contributions and constrained DOFs are removed before solving.

## 4. Generalized eigenproblem

`K φ = λ M φ`,

`ω = sqrt(λ)`, `f = ω/(2π)`.

The sparse solver uses internal sparse matrix storage, RCM reordering when requested, IC(0)-preconditioned conjugate-gradient solves and block inverse/subspace iteration with a dense Rayleigh-Ritz projection. Before iteration, K and M are independently normalized by their maximum diagonal scales; the physical eigenvalue is restored by the K/M scale ratio. This removes dependence on absolute matrix coefficient magnitude without changing the generalized eigenproblem. The M-orthogonalization uses two-pass modified Gram-Schmidt. Eigenpair residuals and positivity checks are validation requirements. Frequency order is not a physical mode identifier near crossings or degeneracies.

## 5. Modal tracking

Three different implementations have distinct roles:

- legacy optimizer shape tracking reorders `ModalData` using the available profile traces;
- cylindrical FEM tracking computes mass-weighted overlaps and subspace principal-angle similarities in a common DOF space, and can exactly embed nested Fourier bases when the meridional DOFs are unchanged;
- the diagnostic state tracker predicts frequency with a constant-velocity Kalman state `[f, df/ds]`, then associates candidate modes using boundary-shape/subspace evidence and frequency likelihoods.

The diagnostic tracker's full-field similarity currently uses row-sum lumped mass weights only when node layouts match. Cross-mesh volume interpolation, physical arc-length registration and consistent-mass principal-angle comparison are not yet implemented. The filter predicts associations; it does not modify FEM eigenpairs.

## 6. Excitation

The current optimizer uses a meridional generalized influence proxy. The Chinese-bell strike model is separate: it constructs a spatially localized nodal force distribution along the local surface normal. For a mode `φ` and force vector `f`, the generalized force is

`Q = φᵀ f`.

A mass-normalized participation measure divides `|Q|` by `sqrt(φᵀ M φ)`. This force projection must not be replaced by `φᵀ M f`.

Strike position and width are physical inputs. Lengths in the direct Chinese strike interface are SI metres; angular widths are radians. Contact width is not an optimizer knob absent measured striker/contact evidence.

## 7. Acoustic transfer and synthesis

The current production synthesis uses a reduced acoustic modal representation. The Chinese objective evaluator includes directional monopole/dipole radiation proxies and a decay gate, but these are not a calibrated non-axisymmetric radiation solution. Mechanical modal amplitude, generalized force, radiated pressure and perceived loudness are different quantities.

## 8. Coordinate morphing status

`src/math/modal_coordinate_morphing.*` remains active in the legacy optimizer and geometry-listening workflow. It is a reduced-order Rayleigh-energy estimate applied to a coordinate-morphed mode, not a new FEM solution on the morphed geometry. The direct cylindrical FEM path is separate and currently diagnostic/opt-in. The morphing code must not be removed until its active call sites are migrated and regression behavior is deliberately replaced.

## 9. Scale law

For uniform geometric scale `s` at fixed elastic modulus and density, `K(s)∝s` and `M(s)∝s³`, hence `f(s)∝1/s`. This provides a useful dimensional regression check; it does not replace physical material calibration.

## 10. Model limitations

The code does not currently claim full 3-D structural FEM, calibrated near-field acoustics, validated impact/contact mechanics, complete cross-mesh modal-field continuation, or historical authenticity of generated geometries. The Chinese objective remains outside production optimization until its mechanical and acoustic gates pass.

## 11. Acceptance status and numerical gates

The project distinguishes three evidence states:

- **VALIDATED** — the metric, threshold, reference and automated verification are defined and exercised by a reproducible test.
- **PROVISIONAL** — the metric is implemented or experimentally useful, but its threshold or reference is not yet validated.
- **EXPERIMENTAL** — a research/diagnostic evaluator that must not be interpreted as a production physical gate.

The normative list of acceptance criteria is maintained in `docs/VALIDATION_THRESHOLDS.md`. A gate is considered passed only when its automated verification reports PASS. Changing a threshold is an engineering change: the threshold rationale and affected tests must be updated together. A failed candidate is not grounds for weakening a gate.

For the implemented generalized eigensolver, the final residual is evaluated in the normalized problem as

`r = ||K_n x - lambda M_n x|| / ((||K_n|| + |lambda| ||M_n||) ||x||)`.

The denominator is never replaced by an eigenvalue floor. The current convergence contract requires `r < max(10*tolerance, 1e-10)` and M-orthogonality error `< 1e-8`; with the normal `tolerance=1e-8` this means a residual gate of `<1e-7`. These are numerical solver contracts, not experimental accuracy claims.

Near-degenerate eigenvalues are treated as subspaces for modal comparison where the implementation provides principal-overlap metrics. Individual eigenvector identity is not assumed to be physically meaningful at a crossing or degeneracy. The current cross-mesh modal comparison remains PROVISIONAL until physical arc-length registration, consistent-mass equivalence and broader real-FEM continuation coverage are closed.

The Chinese objective remains EXPERIMENTAL. In particular, its current interval configuration is not accepted as a historical acoustic gate: the implementation currently uses a configurable cents target, while the historical interval statement in the acoustic canon is being revised and requires an explicit, sourced reference.


## Chinese strike generalized-force model (1.2.102)

For a finite-width strike on the non-circular outer surface, the mechanical excitation is represented in the same real Fourier basis used by the cylindrical FEM. The generalized force coefficient for basis function `B_b(φ)` is

`F_b = ∫_S g(z,φ) n(z,φ) B_b(φ) dS`,

where `g` is the localized contact weighting and `n` is the physical outward surface normal. For `r(φ,z)=α(φ) r_ref(z)`, both `n` and the surface-area Jacobian depend on `α(φ)`. Therefore they must be integrated jointly over azimuth. Evaluating them only at the strike azimuth and multiplying by a scalar Fourier coefficient is not equivalent for an elliptical section.

Modal participation remains `Q_i = φ_i^T F / sqrt(φ_i^T M φ_i)` for the unit-modal-mass convention. This is a mechanical participation diagnostic; it is not a radiation or perceived-loudness measure.

The two-tone branch decision remains separate from the family-continuation metric. A high subspace similarity only proves that the same two-dimensional family is present; it does not prove that the two physical strike conditions select opposite branches.


## Chinese two-tone robustness and interval continuation (1.2.105)

The mechanical diagnostic now separates three questions: (1) whether a two-dimensional modal family can be continued, (2) whether long- and short-axis generalized-force loads select opposite branches, and (3) whether the resulting branch frequencies reach the documented interval region. The tests use the integrated strike load

`F_b = ∫_S g(z,φ) n(z,φ) B_b(φ) dS`

and M-weighted modal/subspace association. A positive branch-selection observation is not sufficient for physical acceptance: the branch split must also be robust to mesh and strike-height changes and eventually agree with independent interval evidence.

The 1.2.105 diagnostic found a controlled geometry family with robust opposite branch selection around q=0.85 on both 15- and 23-node meshes. Increasing ellipticity raises the split, but q=0.70 reaches only about 206–212 cents, below the documented 1.189–1.260 frequency-ratio region. This indicates that meridional profile shape, not merely ellipticity, must be investigated next.

## Chinese modal-family identification (1.2.104)

The Chinese mechanical audit now treats a modal family as a two-dimensional subspace rather than a frequency-rank pair. Candidate seeds are enumerated in an operational near-circular reference, then associated across ellipticity by the minimum-principal-overlap metric. Individual branch overlaps are reported separately because vectors inside a nearly degenerate family can exchange identity even when the family subspace remains stable.

For each associated pair, the strike response is evaluated from the physical generalized force

`Q_i = phi_i^T f / sqrt(phi_i^T M phi_i)`,

where `f` is the integrated surface-normal strike load. Long-axis and short-axis responses are compared only as a diagnostic branch-orientation observation. No numerical split, similarity or contrast from this matrix is promoted to a physical threshold.

## Meridional profile variation and two-tone split (1.2.106)

The diagnostic matrix keeps the non-axisymmetric cross-section parameter q and the surface generalized-force strike model fixed while varying the normalized meridional radius profile r(z)/R. This separates two mechanisms: cross-sectional symmetry breaking controls branch orientation, while meridional geometry changes the stiffness/mass distribution and therefore the frequency separation of the continued family.

For the tested profiles, the maximum observed split was 356.814 cents at q=0.70. The documented ratio interval 1.189–1.260 corresponds to approximately 299.7–400.1 cents. However, the same q=0.70 state did not retain opposite long/short branch selection. A q=0.75 state of the same family retained opposite branch selection but reached only 251.4 cents. Thus interval matching and strike selectivity must remain separate conditions.


## Refined two-tone mechanical robustness (1.2.107)

The current diagnostic candidate is followed from the near-circular reference through q=0.85, 0.80, 0.75, 0.73, 0.72, 0.71 and 0.70. For each continuation step, the associated two-dimensional modal subspace is retained before projecting the surface-integrated strike force onto its two branches.

The refined 31- and 35-node meshes give frequency splits of approximately 305, 328, 351 and 372–373 cents for q=0.73, 0.72, 0.71 and 0.70 respectively. The same states show opposite branch preference for ranges of strike height. This establishes a diagnostic conjunction between frequency separation and directional excitation in the current mechanical model, but it does not establish historical authenticity or a production acceptance threshold.

The earlier q=0.70 result at 23 nodes produced a much larger branch/split anomaly with poor intermediate continuation similarity. It is therefore treated as a mesh/association warning and is not used as evidence. The refined 31/35-node trend is the relevant current diagnostic result.

## Focused profile/mesh refinement (1.2.111)

The 1.2.111 diagnostic preserves the same continuation topology used in the 1.2.107 candidate: the pair is seeded in the near-circular reference and continued through q=0.85, 0.80, 0.75, 0.73, 0.72, 0.71 and 0.70. Only after that continuation are the long-axis and short-axis surface-integrated strike loads projected onto the tracked two-dimensional family. This avoids the earlier methodological error of jumping directly from q=0.9999 to a strongly deformed state.

The focused 31/35-node comparison shows that the frequency split is comparatively mesh-stable (maximum observed difference <2 cents), but branch selection is not invariant. The base profile retains the q=0.70 conjunction on both meshes; `k1 +0.005` loses it on the 35-node q=0.70 case; `k3 +0.005` loses it on the 35-node q=0.73 case but retains it at q=0.70. Therefore modal-family similarity and directional excitation must remain separate diagnostics from frequency-ratio matching.

## Meridional-profile perturbation sensitivity (1.2.109)

The controlled perturbation audit changes one internal normalized radial knot at a time by ±0.005 while preserving the endpoints. The modal family is seeded in the near-circular reference and continued through the same ellipticity path. The audit demonstrates that modal splitting and branch excitation are both geometry-sensitive; therefore frequency-ratio matching without profile and mesh continuation is not sufficient evidence for a Chinese two-tone candidate.

The perturbation screen is diagnostic only. Its subspace-similarity conjunction is not a physical threshold and must not be used as an optimizer gate.

## 1.2.111 contact-footprint sensitivity

The Chinese strike is a normalized distributed generalized force over the mapped outer surface. Meridional and angular Gaussian widths control the spatial weighting of that load and therefore modal participation, but they do not alter the structural eigenproblem. The 1.2.111 diagnostic varies these widths independently while keeping the corrected surface-normal/surface-area integration and the established modal-family continuation path. Because no external striker/contact measurement currently constrains these widths, the scenarios are sensitivity probes rather than physical parameter estimates.


## 1.2.112 contact-footprint parameterization
The localized strike is a distributed generalized-force integral over the outer surface. In the current implementation the meridional Gaussian is expressed in z and the angular Gaussian in phi. For a local physical tangential width `w`, a first-order local mapping is `w_z = w/(ds/dz)` along the meridian and `w_phi = w/r` along the circumferential direction. These transformations are only a parameterization audit, not a claim about the actual striker shape or contact mechanics. The 35-node q=0.70 diagnostic produced relative normalized-influence differences of roughly 0.17–0.19 between the physical-width mapping and the untransformed coordinate-width interpretation. Therefore contact geometry must be established independently before widths can be treated as physical quantities.


## Contact footprint surface metric (1.2.113)
For the mapped surface `R(phi,z)=alpha(phi) r_ref(z)`, a local physical length must use the surface metric rather than raw coordinate widths. With `alpha=1/sqrt(cos^2(phi)+sin^2(phi)/q^2)`, the meridional factor is `g_z = sqrt(1 + alpha^2 (dr_ref/dz)^2)` and the circumferential factor is `g_phi = alpha r_ref sqrt(1 + (d log(alpha)/dphi)^2)`. At the principal axes `d log(alpha)/dphi=0`, so `g_phi=r_ref` on the long axis and `g_phi=q r_ref` on the short axis. This corrects the first-order 1.2.112 interpretation. The corrected mapping is used only in diagnostics; the production strike API remains expressed in z and phi until independent striker/contact evidence justifies an API change.
