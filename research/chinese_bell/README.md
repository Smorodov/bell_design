# Chinese-bell research module

This directory is a **detached research branch** for the Chinese two-tone bell investigation.

It is intentionally excluded from the normal BellOptimizer build. The production project does not compile or link this module unless `BUILD_CHINESE_RESEARCH=ON` is explicitly requested.

## Why it is detached

The investigation established useful numerical machinery, but the current reduced FEM does not reproduce the physically documented Chinese `(m,n)a/(m,n)b` modal-family mechanism with sufficient confidence.

Keeping the code here has two purposes:

1. preserve the work so another researcher can continue without restarting from zero;
2. prevent unresolved experimental assumptions from becoming part of the normal production build.

The correct interpretation is **research infrastructure, not a production feature**.

## What is preserved

- corrected non-axisymmetric surface-normal strike integration;
- finite-contact surface metric handling;
- temporal impulse diagnostic;
- modal-family/subspace continuation experiments;
- modal topology diagnostics;
- Chinese two-tone objective experiments;
- all focused Chinese research contracts and regression tests;
- consolidated scientific conclusions and literature references.

## What is deliberately not claimed

The current module does not establish:

- a physical `(m,n)a/(m,n)b` modal family;
- robust reciprocal `A/B` strike selectivity;
- a production Chinese-bell frequency target;
- a validated Chinese acoustic/radiation/damping objective;
- a physical striker/contact law;
- historical authenticity of the reduced model.

## Build

Normal BellOptimizer builds do **not** compile this directory.

To explicitly build the research module from the project root:

```sh
cmake -S . -B build-chinese -DCMAKE_BUILD_TYPE=Release -DBUILD_CHINESE_RESEARCH=ON
cmake --build build-chinese --parallel
ctest --test-dir build-chinese --output-on-failure
```

The module consumes the normal `bell_core` target. It does not modify that target with Chinese-specific sources.

## Reopening the research

Read these documents first:

1. `docs/CHINESE_RESEARCH_HANDOFF_1.2.126.md`
2. `docs/CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`
3. `docs/CHINESE_MODAL_TOPOLOGY_AUDIT.md`
4. `docs/CHINESE_MODAL_FAMILY_CONTRACT.md`
5. `docs/CHINESE_BELL_LITERATURE.md`

Do not start by tuning tracker thresholds, ellipticity, contact width, or frequency windows. The documented failure modes are precisely the traps this separation is intended to prevent.

## Reusable ideas for the main project

The detached code is not a dead end. Its numerical work can be reviewed independently for reuse in the production circular/elliptic bell path, especially:

- exact surface metrics;
- local surface-normal generalized forces;
- subspace rather than frequency-rank continuation;
- full-field modal comparison.

Such reuse must be re-derived against the production mechanism rather than copied wholesale.
