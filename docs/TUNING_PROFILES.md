# Acoustic target profiles

## Architectural rule

A tuning profile is an **acoustic target data set**, not a separate optimizer implementation. Western minor-third, Western major-third, neutral, the experimental Sysoi reference and user-defined profiles must pass through the same generic objective evaluator. Profile names and cultural labels are descriptive metadata; they must not select hard-coded scoring branches.

Geometry priors are independent. `seed_family` belongs to `evolution*.json`; the acoustic target file must not silently choose a geometry family. See [`TARGET_PROFILE_ARCHITECTURE.md`](TARGET_PROFILE_ARCHITECTURE.md) for the proposed generic JSON schema and migration requirements.

## Existing profile catalogue (legacy JSON v1)

Stage 10 accepts both legacy v1 and generic v2 profiles. The v1 loader preserves the four-slot mapping `{0,2,3,4}` relative to reference mode `1`, with amplitude targets `{1,2,3,4}` relative to mode `0`. V2 allows variable-length targets and explicit sorted-index selectors. These remain sorted-index proxies, not a physical mode-family identification system.

Files currently present under `bells/`:

- `targets_european.json` — Western minor-third reference, 440 Hz.
- `targets_european_major_third.json` — Western major-third variant, 440 Hz.
- `targets_neutral.json` — neutral geometry label with the legacy minor-third acoustic ratios, 440 Hz.
- `targets_western_220.json`, `targets_western_440.json`, `targets_western_880.json` — Western minor-third reference at three reference frequencies.
- `targets_orthodox.json`, `targets_orthodox_220.json`, `targets_orthodox_440.json`, `targets_orthodox_880.json` — experimental Sysoi single-exemplar ratios at the indicated F2 reference frequency.
- `targets.json` — default legacy Western reference.

The corresponding `evolution*.json` files select the geometry seed family and search parameters. `evolution.seed_family` takes precedence; the legacy target `tradition` key is only a fallback when no seed family is specified. Conflicting values are rejected.

## Western reference

The minor-third reference is:

`Hum : Prime : Tierce : Quint : Nominal = 0.5 : 1.0 : 1.2 : 1.5 : 2.0`

The major-third variant is:

`Hum : Prime : Tierce : Quint : Nominal = 0.5 : 1.0 : 1.25 : 1.5 : 2.0`

These are reference profiles, not a claim that all Western bells share one exact spectrum.

## Russian reference — experimental Sysoi exemplar

The experimental target is based on Table 1 for the Sysoi bell at Rostov Uspensky Cathedral in I. Aldoshina and A. Nicanorov, “The Investigation of Acoustical Characteristic of Russian Bells,” AES 108th Convention preprint 5117 (2000):

- F1 = 58.5 Hz
- F2 = 115.5 Hz (reference)
- F3 = 145 Hz
- F4 = 189 Hz
- F5 = 252.5 Hz

Ratios relative to F2 are approximately:

`F1/F2 : F3/F2 : F4/F2 : F5/F2 = 0.51 : 1.26 : 1.64 : 2.19`

This is one measured historical exemplar, not a universal Russian Orthodox tuning law. The paper reports that familiar Western octave/minor-third/fifth relations were not identified in the examined Russian bells and discusses differences across bell size groups. More measured Russian bells are required before generalizing the target.

The current executable treats sorted FEM slots as provisional F1–F5 proxies. Their physical correspondence to measured spectral components is not independently established. Amplitude objective weights are zero in the experimental JSON because published SPL values depend on strike and observation time and are not transferable modal-amplitude targets. The profile remains exploratory and does not certify historical or liturgical authenticity.

## User-defined profiles

V2 supports variable-length `frequency_targets` and `amplitude_targets`, an explicit `frequency_reference`, and explicit amplitude reference modes. A working example is `examples/target_profiles/custom_v2_example.json`. Invoke it with `bell optimize --targets examples/target_profiles/custom_v2_example.json --evolution bells/evolution.json`.

Only `sorted_index` selectors are implemented. Indices are zero-based sorted modal slots; they do not establish physical mode-family identity, especially for near-degenerate circular-bell pairs. `profile_id` and `display_name` are descriptive metadata and do not select a scoring branch. The v1 files remain supported through a conversion layer. Full numerical ranking parity with the pre-refactor Western optimizer is still pending.

## Chinese zhong

Chinese two-tone targets remain in the separate deferred research branch. The current production axisymmetric model does not represent the required non-axisymmetric modal families; the Chinese research JSON is not a supported production tuning profile.
