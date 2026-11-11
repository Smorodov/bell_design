# Development and release workflow

## Checkpoint discipline

- Start from the newest archive actually present. Verify filename, ZIP integrity and project version before editing.
- Keep the input archive unchanged; work in a separate extracted directory.
- Do not claim an artifact exists until its exact path has been created and checked.

## Validation ladder

1. Inspect the affected source, call sites and tests.
2. Build the changed production target and directly affected contract tests.
3. Run focused CTest contracts.
4. Run full build/CTest for cross-cutting changes or release candidates.
5. Run mesh/mode/time-step convergence or physical calibration only when the engineering decision requires it.

Record the actual command, configuration, test result and any unrun checks. A focused pass is not a full-suite pass.

## Research and data retention

Keep source probes, tests and the minimum reproducibility data needed to support durable conclusions. Replace transient WAV bundles, duplicate CSV sweeps and raw exploratory dumps with compact tables and documented methods when the raw files are no longer needed. Keep Chinese-bell research isolated from the active round-bell production path.

## Code-change discipline

Prefer one coherent engineering increment. Do not combine physical-model changes with unrelated refactoring. Never weaken numerical tolerances or normalize physical intermediate quantities merely to make a preview sound better. Listening renders are diagnostics unless the entire acoustic chain is calibrated.

## Release checklist

- Confirm project version and release note.
- Run the appropriate build/tests.
- Verify documentation links and remove stale stage-specific documents.
- Exclude build directories, caches, temporary logs, WAVs and generated test output from the archive.
- Inspect archive contents and run ZIP integrity verification.
- Record SHA-256 and provide only the verified artifact path.
