# Modal analysis and tracking

## Three separate tasks

- Legacy optimizer modal ordering (`core/modal_tracking.*`) uses profile-shape information available in `ModalData`.
- Direct vector-FEM same-space comparison (`fem/cylindrical_modal_tracking.*`) uses mass-weighted overlaps, principal angles and exact nested Fourier-basis embedding.
- Diagnostic continuation (`core/modal_state_tracker.*`) uses a constant-velocity Kalman frequency predictor plus shape/subspace evidence. It is test-only and not wired into the optimizer.

Cross-mesh full-Q4 field comparison (`fem/cylindrical_modal_field.*`) remains a separate diagnostic utility. It integrates each Q4 field on its own mesh over the common physical domain. The state tracker does not own mesh topology, so merging these APIs would obscure their distinct responsibilities.

## Completed continuation evidence

- Same-geometry q sequence 0.99→0.96: 4/4 tracker updates accepted; minimum pair similarity 0.999963; no ambiguous updates.
- 15→23 profile-node continuation: tracker accepted the update with pair similarity 0.999909. The legacy boundary-shape matcher did not report shape tracking on that case; this is not evidence that the tracker itself uses the full cross-mesh Q4 metric.
- Fourier enrichment `{0,2,4}`→`{0,2,4,6}`: pair similarity 0.999960.
- Controlled 1.2.159 cross-mesh assignment patch: minimum individual overlap 0.999854; 2D principal subspace similarity 0.999966; one-to-one assignment preserved all four reference identities. This is a controlled regression, not broad real-bell generalization.

## Acoustic mesh convergence evidence

For one q=0.85 real-bell case and tracked two-mode family 10:11, profile nodes 15/23/31/41/49 gave:

| Refinement | Max pair frequency change vs 15 nodes | Pair radiation norm change vs 15 nodes |
|---|---:|---:|
| 23 | 2.2256% | 4.0712% |
| 31 | 2.9786% | 5.4739% |
| 41 | 3.3691% | 6.2112% |
| 49 | 3.5233% | 6.5045% |

The adjacent 41→49 refinement changed maximum pair frequency by 0.1596% and pair radiation norm by 0.3127% at each ear. The coarse 15-node result is therefore not accepted as acoustically converged for this metric. Meridional quadrature 8→16 on the fixed profile changed the complex field by 0.241344%.

## Closure and limitations

The diagnostic tracking stage is closed for the tested q path, selected mesh/Fourier continuations and selected modal families. This is not a universal proof for arbitrary bells, crossings or spectra. Particle filters, multi-hypothesis tracking and a large-spectrum assignment API are deferred until a concrete case demonstrates that the current method fails.

Frequency rank alone is never a physical modal identity criterion. Subspace continuity is a numerical association measure, not proof of a historical mode family. No tracker is integrated into the production optimizer by this audit.
