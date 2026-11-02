# Chinese Two-Tone Bell — Consolidated Physical and Numerical Audit

## Purpose

This document replaces the sequence of exploratory Chinese-bell audit notes accumulated during 1.2.102–1.2.121.

The historical sweep data, individual parameter matrices, intermediate logs, and experiment-specific narratives are intentionally not retained here. Their conclusions are consolidated below so that the project has one durable engineering record rather than a chain of experimental modifications.

The Chinese-bell implementation remains **diagnostic / experimental**. It is not a production optimizer objective and none of the conclusions below constitute a validated historical-authenticity claim.

## 1. Physical facts established by external evidence

The literature consistently supports the following points:

1. Ancient Chinese two-tone bells can produce two distinct tones from one bell at different strike positions.
2. The bells have non-circular / almond- or lens-shaped cross sections.
3. Vibrational modes occur in pairs with closely related spatial structure; paired members differ in the placement of nodes and antinodes.
4. Different strike positions selectively excite the two members of the relevant modal pair.
5. The interval is commonly associated with a minor or major third, but the project does not treat a broad numerical interval band as sufficient evidence of the physical mechanism.
6. Chinese bells also exhibit comparatively rapid attenuation, which is acoustically relevant but is not yet a validated production target in this project.

Primary sources are listed in `CHINESE_BELL_LITERATURE.md`.

## 2. What the project model establishes numerically

### 2.1 Structural FEM core

The generalized eigenproblem and sparse structural pipeline are numerically established and covered by the general BellOptimizer validation suite. The Chinese work does not alter the structural solver architecture.

### 2.2 Non-axisymmetric strike excitation

The strike model was corrected to integrate the physical surface normal over the localized load instead of multiplying a scalar Fourier coefficient by a single normal vector.

The current formulation therefore uses the local elliptical surface geometry and its surface metric during generalized-force assembly.

This is a numerical/modeling correction, not an experimental claim about the real striker.

### 2.3 Contact-width semantics

The model parameters formerly explored as 20/35/50 mm are **mathematical spatial-width scenarios** only.

They are not identified as:

- mallet diameter;
- contact-patch diameter;
- FWHM;
- Gaussian standard deviation of a measured physical pressure field.

Available literature does not provide enough striker/contact measurements to justify any of these identifications for the present model.

### 2.4 Temporal impact

A normalized temporal impulse layer has been tested diagnostically. It demonstrates the expected mathematical relationship between shorter impact duration and broader spectral excitation.

No universal impact duration, velocity, contact stiffness, Hertz exponent, or force-time history has been promoted to production parameters.

## 3. Modal-family investigation — consolidated result

The project investigated several near-degenerate modal families from the nearly circular reference configuration and followed them through controlled ellipticity changes using mass-weighted two-dimensional subspace continuation.

The important numerical lesson is that **frequency rank is not a reliable modal-family identifier**.

An early narrow modal window produced apparent 300–370 cent separations. The wider 16-mode audit showed that the apparent large split was associated with a branch/window-selection artifact rather than a robust continuation of the original reference pair.

The corrected continuation of the selected reference family remains highly similar as a two-dimensional subspace, but its frequency separation stays comparatively small over the tested ellipticity range.

Therefore:

> A large frequency split obtained after modal re-association is not accepted as evidence for the Chinese two-tone mechanism.

## 4. Spatial modal topology — consolidated result

The spatial structure of the tracked pair was then examined using surface normal displacement rather than frequency alone.

The physical literature indicates that the relevant paired modes should have the same nodal-family orders while their node/antinode placement is complementary, particularly at the spine/xian and the corresponding strike regions.

The current reduced FEM model does not yet reproduce that physical identification robustly for the tracked candidate pair.

Thus the current result is:

- subspace continuation: numerically stable;
- physical `(m,n)a/(m,n)b` identification: **not established**;
- A/B strike complementarity: **not established**;
- historical two-tone mechanism: **not validated**.

## 5. Strike selectivity — consolidated result

Strike projections and local surface-normal modal response were evaluated at the two principal strike directions.

The calculations demonstrate that strike location materially changes modal excitation, but they do not yet produce a robust reciprocal pattern in which:

- strike A strongly excites mode A and suppresses mode B; and
- strike B strongly excites mode B and suppresses mode A.

Consequently, frequency separation alone must not be used as a proxy for strike selectivity.

## 6. Meridional profile and mesh studies

Meridional profile changes materially affect modal splitting and branch selection. Several promising frequency-ratio observations were found, but the same configurations did not simultaneously preserve the required modal-family identity and strike selectivity.

Mesh refinement removed at least one apparently strong split that was later identified as a continuation/resolution artifact.

The durable conclusion is therefore:

> Profile shape, mesh resolution, modal association, and strike selectivity must be treated as coupled validation dimensions.

No profile-specific frequency ratio is promoted to a production target.

## 7. Contact-model studies

The contact audits established three durable conclusions:

1. The exact local surface metric matters for converting a physical-width hypothesis into the model's meridional/angular coordinates.
2. Finite-contact shape can materially change modal generalized forces; therefore a Gaussian footprint is not physically interchangeable with an arbitrary finite patch.
3. There is currently insufficient independent evidence to select a unique striker geometry, pressure distribution, or contact law.

Therefore the project does **not** add a physical striker parameter to the production optimizer at this stage.

## 8. Current evidence status

| Topic | Status | Durable conclusion |
|---|---|---|
| Structural FEM | VALIDATED numerically | Solver and matrix contracts remain intact. |
| Non-axisymmetric surface-normal excitation | VALIDATED numerically | Correct local surface geometry is used. |
| Modal continuation machinery | VALIDATED numerically / physical thresholds provisional | Subspace continuation is preferred over frequency rank. |
| Chinese modal family identity | PROVISIONAL | Current candidate family is not physically identified with sufficient confidence. |
| A/B strike selectivity | PROVISIONAL | Not yet demonstrated robustly. |
| Chinese frequency interval | PROVISIONAL | No production interval gate. |
| Striker/contact geometry | UNDERDETERMINED | Requires independent physical evidence. |
| Temporal impact model | EXPERIMENTAL | Mathematical diagnostic only. |
| Radiation/damping | EXPERIMENTAL | Not yet a validated Chinese production model. |
| Historical authenticity | NOT ESTABLISHED | No heuristic is sufficient for this claim. |

## 9. Decisions carried forward

- Keep Chinese optimization isolated from the Western production objective.
- Do not tune tracker thresholds to manufacture a preferred branch.
- Do not select a modal family by frequency rank alone.
- Do not interpret 20/35/50 mm as a physical striker/contact dimension.
- Do not add Hertz/contact parameters without independent striker/contact evidence.
- Do not promote a minor/major-third interval into a production gate by itself.
- Do not use acoustic proxies to compensate for unresolved mechanical modal identity.
- Do not claim historical authenticity from the current reduced model.
- Keep numerical contract tests; remove experiment-specific sweep executables and raw result documents.

## 10. Remaining scientific question

The next meaningful task is no longer another parameter sweep. It is to determine whether the current reduced FEM representation contains a physically appropriate modal family whose spatial topology and strike selectivity agree with the documented Chinese `(m,n)a/(m,n)b` mechanism.

If that cannot be demonstrated, the project should explicitly record the limitation of the reduced model rather than compensate with additional tuning parameters.
