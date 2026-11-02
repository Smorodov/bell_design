# Chinese Bell Model — Decision and Research Freeze (1.2.124)

## Purpose

This document converts the consolidated Chinese-bell audit into an architectural decision. It is the current decision point for the Chinese-bell research path and replaces the need to revisit the historical experiment sequence.

## Decision

The current reduced Chinese-bell model is **not promoted to a production optimization objective**.

No further parameter sweep of meridional profile, ellipticity, contact width, temporal impact, or modal-tracker thresholds is authorized as a substitute for the missing physical identification.

The current implementation remains diagnostic/experimental and isolated from the Western production objective.

## What is established

- The structural FEM core is numerically validated by the general solver contracts.
- The non-axisymmetric strike projection now uses the local elliptical surface geometry and surface metric.
- Strike position materially changes modal excitation.
- Modal continuation should use mass-weighted subspace identity rather than frequency rank.
- Contact-width values used in diagnostics are mathematical spatial scenarios, not identified physical striker dimensions.
- A temporal impulse layer can be used as a mathematical diagnostic, but its physical parameters are not identified.

## What is not established

The current reduced model has not demonstrated all of the following simultaneously:

1. a physically identified `(m,n)a/(m,n)b` modal family;
2. the documented complementary spatial topology of that pair;
3. reciprocal A/B strike selectivity;
4. a stable frequency interval associated with that same physical pair.

A large frequency separation obtained after modal reassociation is therefore not accepted as evidence of the Chinese two-tone mechanism.

## Literature-based physical contract

The maintained literature supports the following physical picture:

- Chinese two-tone bells can produce distinct tones at different strike positions.
- Relevant modes occur in pairs with closely related nodal structure and complementary node/antinode placement.
- Strike position selectively excites different members of a modal pair.
- The bell cross-section is non-circular and is central to the modal splitting mechanism.

These statements are literature-derived. They are not claims that the current reduced FEM model reproduces them.

## Research gate for reopening the model

The Chinese model should only be reopened for a new numerical development when there is a clearly stated physical question that can distinguish among these possibilities:

- the current reduced shell/FEM representation is insufficient;
- the selected modal family is incorrect;
- the strike representation is insufficient;
- the available geometry is insufficient;
- or the required experimental reference data are missing.

The next useful evidence should therefore be **independent modal-shape/geometry data**, not another unconstrained parameter sweep.

## Explicitly deferred

The following are deferred until the physical modal-family problem is resolved:

- radiation model;
- damping/T60 model;
- acoustic roughness/beat objective;
- production Chinese optimization objective;
- physical striker/contact parameters;
- historical-authenticity claims.

## Cleanup policy

Historical sweep data, experiment-specific executables and raw logs are not part of the maintained project surface. They may remain in external research archives, but they are not required to understand or build the current project.

The maintained documentation consists of:

- `CHINESE_BELL_ACOUSTIC_CANON.md`
- `CHINESE_BELL_STATUS.md`
- `CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`
- `CHINESE_BELL_LITERATURE.md`
- this decision document

## Status

**Chinese two-tone mechanism: RESEARCH BLOCKED pending physical modal-family identification.**

This is not a failure of the structural FEM solver. It is a limitation of the current evidence/model correspondence.


## 1.2.124 modal-family audit

A dedicated diagnostic contract now measures the direct FEM modal subspace, dominant Fourier order, and explicit long/short strike projections. On the canonical reduced model, the tracked candidate pair retained high numerical subspace similarity (~0.928), but its two members did not retain a common dominant Fourier order (m=1 and m=4 in the diagnostic run). The meridional nodal-circle count remains unresolved. Therefore the physical Chinese modal family remains **RESEARCH_BLOCKED**; no production gate or objective was enabled.
