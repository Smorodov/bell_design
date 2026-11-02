# Validation thresholds and acceptance criteria

## Evidence states

- **VALIDATED** — metric, threshold, reference and automated verification are defined.
- **PROVISIONAL** — implementation exists, but the physical/reference threshold is not established.
- **EXPERIMENTAL** — diagnostic/research only; not a production gate.

## Numerical solver gates

### ACT-EIG-001 — generalized eigen residual

`||K_n x-lambda M_n x|| / ((||K_n||+|lambda| ||M_n||)||x||) < max(10*tolerance,1e-10)`.

State: **VALIDATED**.

### ACT-EIG-002 — M-orthogonality

Pairwise M-orthogonality is required below `1e-8` for the validated eigenpair contract.

State: **VALIDATED**.

### ACT-MAT-001 — sparse symmetry

Relative sparse symmetry mismatch must be `<= 1e-10`.

State: **VALIDATED**.

### ACT-MOD-001 — principal overlap range

Principal-overlap values must remain in `[0,1]`.

State: **VALIDATED**.

## Modal continuation gates

The following remain provisional because physical/reference thresholds are not fixed:

- ACT-MOD-010 — cross-mesh frequency convergence;
- ACT-MOD-011 — cross-mesh modal identity / minimum principal overlap;
- ACT-MOD-012 — association margin;
- ACT-MOD-013 — branch retention;
- ACT-MOD-014 — ambiguity handling.

Frequency rank is never a physical acceptance criterion for modal identity.

## Mechanical Chinese-bell gates

The following remain **PROVISIONAL / TBD**:

- ACT-CN-MEC-001 — two-tone modal-family existence;
- ACT-CN-MEC-002 — strike selectivity;
- ACT-CN-MEC-003 — pair stability;
- ACT-CN-MEC-004 — contact-model domain;
- ACT-CN-MEC-005 — contact-width semantic identity.

No numerical threshold is promoted to a physical gate until an independent reference or measurement is identified.

The consolidated evidence and negative results are maintained in `research/chinese_bell/docs/CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`.

## Acoustic gates

| ID | Gate | Metric | Threshold | State |
|---|---|---|---|---|
| ACT-AC-001 | Radiation model | agreement with independent reference | TBD | EXPERIMENTAL |
| ACT-AC-002 | Perceived pitch | named psychoacoustic metric | TBD | EXPERIMENTAL |
| ACT-AC-003 | Roughness/beating | named psychoacoustic metric | TBD | EXPERIMENTAL |
| ACT-AC-004 | Decay | T60 or band-limited decay | TBD | EXPERIMENTAL |
| ACT-AC-005 | Pair separation | joint separation/competition metric | TBD | EXPERIMENTAL |

Current monopole/dipole proxies are diagnostic only.

## Chinese interval gate

The minor/major-third interval is an externally documented characteristic of Chinese two-tone bells, but the project does not use a broad numerical interval band as a sufficient mechanical gate. No production interval threshold is declared until the physical modal family and strike selectivity are established.

## Candidate disposition

A candidate is accepted only when every required gate for its optimization profile passes. Mechanical PASS with acoustic FAIL is REJECTED with explicit failing gate IDs. Experimental metrics cannot silently become PASS criteria.

## CalculiX reassessment trigger

CalculiX remains a differential reference while the declared internal-vs-external validation matrix is incomplete. Reassessment is triggered only after complete declared coverage exists for the current geometry families, ellipticity cases and modal families.

## Morphing-to-FEM migration

Removal of legacy coordinate morphing requires active-call-site migration, direct-FEM regression coverage, defined equivalence thresholds and documentation of intentional differences.

## Non-negotiable rules

- No tuning of gates to make a candidate pass.
- No weakening of a gate solely because an optimization run fails it.
- No silent mixing of Western and Chinese objectives.
- No experimental evaluator as a production physical gate.
- No historical-authenticity claim from geometric or acoustic heuristics alone.
