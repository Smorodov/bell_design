# Development workflow and resource budget

This document defines how development sessions are organized so that limited analysis, tool and compute budgets are spent on new evidence rather than repeated work. It applies to the project and to collaboration sessions that modify or validate it.

## 1. Start from the newest verified checkpoint

- Use the newest archive that is actually available in the current conversation or Library. Confirm its filename, version, size and archive integrity before editing.
- Do not assume that a version mentioned in chat exists as a downloadable artifact. A release exists only after its file has been created and its exact path verified.
- Keep the input archive unchanged. Work in a separate directory and produce a new versioned archive.
- Before editing, inspect the relevant source, tests, documentation and release notes. Avoid broad scans when a focused search answers the question.

## 2. Plan work as small, reviewable increments

For each increment, write down: (a) the question or defect, (b) the smallest change that can answer it, (c) the evidence required, and (d) the stop condition. Prefer one coherent change per release. Do not mix a physical-model change, refactor and unrelated documentation rewrite unless they are inseparable.

When several independent checks can be performed on the same build or dataset, batch them into one planned pass. Do not batch unrelated risky changes merely to reduce the number of sessions.

## 3. Use a validation ladder

Run only the level justified by the change, then escalate when the evidence or risk requires it:

1. **Static inspection / focused diff** — check interfaces, units, assumptions, ownership/lifecycle and affected call sites. No build is needed for a documentation-only change.
2. **Focused contract tests** — run tests directly covering the changed mechanism.
3. **Targeted regression set** — include adjacent mechanisms and known sensitive cases.
4. **Full build and full test suite** — use for cross-cutting changes, release candidates, solver changes or when lower-level results indicate wider risk.
5. **Expensive convergence / physical validation** — run only when the change affects numerical results or when the next decision depends on new convergence evidence.

A test result may be reused only if the relevant source, build configuration, compiler/toolchain, test data and runtime conditions have not changed. Record what was reused and why. If any relevant input changed, rerun the affected level.

## 4. Avoid redundant calculations

- Inspect existing logs, metadata, CSV sidecars, WAV files and test reports before regenerating them.
- Reuse cached intermediate data only when its provenance is known: source revision, parameters, mesh, solver settings and output format. Unknown-provenance data is exploratory, not validation evidence.
- Prefer a cheap diagnostic run to choose the next expensive run; do not treat a coarse run as proof of fine-mesh convergence.
- Use parameter sweeps only to answer a stated sensitivity or convergence question. Stop when the decision criterion is met or when the trend shows that the current model cannot answer the question.
- Do not rerun a successful unchanged test suite merely to produce another identical log.

## 5. Respect tool and analysis quotas

- Tool quotas and reset schedules can vary by plan and product; do not assume a fixed public number of analysis runs or a guaranteed reset time. Use any limit/reset information shown by the product UI as authoritative.
- Prepare a compact analysis batch before invoking a limited tool: define the inputs, requested measurements, outputs and decision that the results will support. Avoid exploratory calls with no decision attached.
- Prefer one well-scoped request that returns all needed measurements over multiple near-duplicate requests. Keep source files and relevant context together to avoid repeated uploads and re-explanation.
- If a required tool reports that its quota is exhausted or unavailable, stop that tool-dependent task. Do not attempt to bypass the limit through another execution route. Continue only with independent work that does not require the unavailable capability, or resume when access returns.
- Clearly distinguish completed analysis, partial analysis, estimates and unrun work. Never infer that a test or measurement passed because it could not be run.

## 6. Preserve reproducibility and evidence

For every meaningful numerical run, retain a concise record of:

- project version / source revision and build configuration;
- exact command, parameters, units and input-data identifiers;
- mesh size, mode count, time step, solver tolerances and relevant environment;
- runtime and exit status;
- key metrics, acceptance thresholds and pass/fail/blocked status;
- output paths and whether the result is diagnostic, converged, calibrated or externally validated.

Keep raw data only when it is needed to reproduce or interpret a result. Summarize durable conclusions in the authoritative documentation; move bulky exploratory or failed sweeps to a clearly isolated research area, and do not include them in normal builds. Never delete the only copy of evidence needed to support a conclusion.

## 7. Release and archive discipline

- Do not create a new downloadable archive after every tiny edit. Accumulate changes into a coherent increment with a clear acceptance condition.
- Before packaging, check version consistency, documentation links, generated-file policy, archive contents and absence of temporary build/cache files.
- Run the validation level appropriate to the increment. For documentation-only changes, perform link/name consistency checks and archive integrity checks; do not claim code tests were run.
- Include a release note stating what changed, what was tested, what was not tested and what remains open.
- Create the archive only after the checks finish; verify that the exact output file exists and can be listed/opened. Provide that exact artifact, not a guessed link.

## 8. Session handoff

At the end of a work session, leave a short handoff containing:

1. current verified project version and archive name;
2. changes completed;
3. commands/tests actually run and their results;
4. checks skipped or blocked, with reasons;
5. next single recommended task and its acceptance condition.

The handoff should be sufficient to resume without rereading the entire conversation. Durable project decisions belong in project documentation, not only in chat.

## 9. Physical-model integrity takes priority over speed

Resource savings must not weaken physical correctness. Do not replace missing evidence with assumptions, relax tolerances solely to obtain a pass, suppress natural beating, trim transients, add filtering or normalize intermediate physical quantities to make a preview sound more pleasing. Listening previews are diagnostic outputs, not calibrated sound-pressure predictions. Preserve raw diagnostic audio and metadata unless a separate, justified presentation-rendering stage is introduced.

For the current internal-clapper work, keep the physical simulation, modal state, acoustic rendering and optimizer integration as separate gates. Do not switch the production optimizer until dimensional scale handling, measured clapper/pivot inputs, modal and audio convergence, contact-model limits and objective acceptance criteria have been addressed.
