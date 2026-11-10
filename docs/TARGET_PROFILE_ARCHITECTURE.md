# Acoustic target profiles — architecture and migration contract

**Status:** implemented baseline in Stage 10. The CLI accepts legacy v1 and selector-based v2 profiles. The first v2 implementation supports only `sorted_index`; modal-family and degenerate-subspace selectors remain future work. Full optimizer ranking parity against a frozen pre-Stage-10 baseline is still an acceptance gate.

## 1. Responsibilities must be independent

A bell run combines four independent inputs:

1. **Acoustic target profile** (`targets*.json`): the desired modal frequency and, when justified, amplitude relationships.
2. **Geometry/search configuration** (`evolution*.json`): seed family, population, generations, mutation, deterministic seed, geometry constraints and material preset.
3. **Synthesis configuration** (`synthesis.json`): sample rate, duration, decay, normalization, receiver and synthesis options.
4. **Solver/model capability**: what modal families and physical quantities the selected solver can actually represent.

A cultural/traditional label is metadata and/or a geometry prior. It must not select a different hard-coded acoustic objective. A profile called `russian_orthodox_sysoi_2000_experimental` must be scored by the same generic target evaluator as a Western, neutral or user-created profile.

The geometry seed family belongs in `evolution*.json` (`seed_family`), not in the acoustic target profile. During migration, the legacy `tradition` key in target JSON may be read as a fallback only when `seed_family` is absent; conflicting values must be rejected or reported clearly, never silently resolved.

## 2. Intended generic JSON shape

The following is an accepted v2 CLI format example:

```json
{
  "schema_version": 2,
  "profile_id": "custom_440_example",
  "display_name": "Custom five-mode example",
  "source": {
    "kind": "user_defined",
    "description": "Illustrative target, not a measured bell"
  },
  "frequency_reference": {
    "mode": { "selector": "sorted_index", "index": 1 },
    "target_hz": 440.0
  },
  "frequency_targets": [
    { "mode": { "selector": "sorted_index", "index": 0 }, "ratio": 0.51, "weight": 1.0, "label": "F1" },
    { "mode": { "selector": "sorted_index", "index": 2 }, "ratio": 1.26, "weight": 1.0, "label": "F3" },
    { "mode": { "selector": "sorted_index", "index": 3 }, "ratio": 1.64, "weight": 1.0, "label": "F4" },
    { "mode": { "selector": "sorted_index", "index": 4 }, "ratio": 2.19, "weight": 1.0, "label": "F5" }
  ],
  "amplitude_targets": []
}
```

### Index and ratio semantics

- Mode indices in this schema are zero-based, refer to the sorted modal list, and are explicitly provisional selectors. They are not claims of physical mode-family identity.
- `frequency_reference.mode` selects the mode whose physical frequency is scaled to `target_hz`.
- Each frequency target is scored as `f(selected_mode) / f(reference_mode)` against `ratio`.
- An amplitude target must explicitly name both its measured mode and its reference mode. Its ratio is the amplitude of the selected mode divided by the reference mode amplitude. No hidden assumption that all amplitudes are relative to the first or “Hum” mode is allowed.
- `weight` must be finite and nonnegative. A zero weight disables that term without deleting its target record.
- Empty target arrays are valid only if at least one active objective term remains. The current evaluator checks selected indices against the available modal count at candidate evaluation time.
- Duplicate selectors, reference/target self-comparisons, invalid indices, non-positive ratios, negative/non-finite weights, and missing required fields are configuration errors.

## 3. Selector evolution

The first implementation may support `sorted_index` only, because that is the selector the current production optimizer can actually evaluate. The schema reserves a selector abstraction so that later validated selectors can include tracked modal families or degenerate modal groups. Do not implement a selector name that the solver cannot honor.

For circular bells, near-degenerate pairs require subspace/group-aware treatment. A profile must not pretend that an arbitrary eigenvector inside a degenerate subspace is a stable named partial. Group selectors are a later capability gate, not something inferred from a JSON label.

## 4. Backward compatibility

Legacy schema v1 remains accepted:

- `prime_freq_target`;
- four-element `freq_ratios` and `freq_weights` arrays;
- four-element `amp_ratios` and `amp_weights` arrays;
- optional `tuning_model` and `tradition` strings.

The v1 compatibility reader translates frequency targets to sorted indices `{0,2,3,4}` relative to index `1`; amplitude targets map `{1,2,3,4}` relative to index `0`. This preserves the historical slot mapping while converting to generic in-memory vectors.

Migration requirements:

1. Keep a v1 reader and translate v1 values into the generic in-memory target representation.
2. Preserve numerical behavior of the production Western baseline exactly within a regression tolerance before enabling v2 by default.
3. Keep v2 parsing strict: wrong types, missing fields, unsupported selectors, non-finite/non-positive ratios, negative weights and duplicate selectors must fail. The first implementation reports the field or target category; deeply nested path diagnostics can be improved separately.
4. Treat `tuning_model` as deprecated descriptive metadata, not an allow-list or objective switch. Stage 10 removes its allow-list gate and profile-specific scoring/reporting branches from the production CLI and optimizer.
5. Keep `seed_family` in the evolution configuration. Read legacy `tradition` only as a documented migration fallback.
6. Do not silently reinterpret old profiles or change the production default objective.

## 5. Profile catalog policy

- Western minor-third, Western major-third and neutral references remain independent profile files.
- The Sysoi-derived Russian reference remains a clearly labelled single-exemplar experimental profile. Its ratios are approximately `0.51, 1.26, 1.64, 2.19` relative to F2, based on the cited published measurements. It is not a universal Orthodox tuning law.
- A profile's cultural name does not guarantee the current axisymmetric solver can represent every physical mode in the measured spectrum.
- Frequency targets and geometry priors are independent: a user may combine a custom acoustic profile with a neutral or other supported geometry seed family.
- Chinese two-tone research remains isolated and deferred; it must not be represented as a supported production target profile until the necessary non-axisymmetric modal families are supported and validated.

## 6. Required implementation and acceptance tests

Before declaring generic profiles implemented:

1. Unit-test v1-to-generic conversion and compare legacy Western objective values against the pre-migration baseline.
2. Load v2 profiles with 1, 4, 5 and more frequency targets; verify every target contributes exactly once.
3. Test arbitrary reference mode indices and amplitude reference selectors.
4. Reject malformed, truncated, overlong, duplicate and out-of-range selectors/weights/ratios with clear diagnostics.
5. Verify a profile's `profile_id` or `display_name` cannot change the objective; only its target data can.
6. Verify changing `evolution.seed_family` changes geometry initialization without changing acoustic targets, and changing the target profile does not silently change the geometry seed family.
7. Unsupported selector capabilities fail during loading; insufficient modal count fails the candidate with the invalid-objective sentinel. A dedicated mode-count regression remains to be added.
8. Test circular near-degenerate mode handling separately; do not claim that sorted-index targets solve modal-family identity.
9. Compare the old and new production Western run using the same geometry seed, deterministic random seed, solver, mesh and evolution settings. Stage 10 preserves the old formula and mapping in the v1 conversion, but full run/ranking parity has not yet been demonstrated. Keep the Western profile as the production reference until this gate passes.

## 7. Documentation rule

Every user-facing guide must distinguish:

- what the current executable accepts today (legacy v1 and selector-based v2 JSON);
- what v2 currently supports (`sorted_index`) versus future selector types;
- what modal selectors the current solver can validate;
- what is only a research/example profile.

Never describe arbitrary-length profiles as implemented until the parser, evaluator, CLI validation and regression suite all support them.
