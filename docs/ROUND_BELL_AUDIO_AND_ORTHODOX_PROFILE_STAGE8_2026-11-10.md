# Round-bell audio demos and experimental Russian profile — Stage 8

Date: 2026-11-10  
Project version: 1.2.184 (unchanged)

## Scope

- Each of the six Western/Orthodox × 220/440/880 Hz launchers now optimizes a case and, on success, creates two WAV previews in that case's `audio/` directory.
- Western targets and the production `western_five_partial` objective are unchanged.
- Orthodox targets are no longer copies of Western ratios. They use a clearly labelled experimental single-bell reference from a published measurement.
- Chinese-bell research remains out of scope.

## Audio outputs

`run_bell_case.cmd` performs the optimization first and requires `best_geometry.json`. It then runs:

1. `bell synth-geometry ...` → `*_modal.wav`. This is a modal listening preview based on a fresh internal modal solve from the optimized geometry. It does not run the optimizer again. The command's existing 15-node mesh and simplified modal excitation/amplitude assumptions remain in force.
2. `bell synth-contact ...` → `*_contact.wav`. This is a localized-strike preview with contact integration and directional radiation, using a 31-node mesh for the demonstration runtime.

The files represent different excitation models and are comparative previews, not calibrated SPL predictions. Contact metadata is written beside the WAV; the modal command also writes envelope/spectrum and morphing metadata files according to its existing behavior.

The internal optimizer may not create `best_bell.dat`. That file is copied only when the solver provides a source modal-data file. The previous guide text implied it was guaranteed; this was corrected.

## Experimental Orthodox target

Reference: I. Aldoshina and A. Nicanorov, “The Investigation of Acoustical Characteristic of Russian Bells,” AES 108th Convention, Paris, 19–22 February 2000, preprint 5117.

Table 1 identifies the Sysoi bell at Rostov Uspensky Cathedral. It lists F1=58.5 Hz, F2=115.5 Hz, F3=145 Hz, F4=189 Hz, F5=252.5 Hz. With F2 used as the prime reference, the target ratios are:

- Hum slot F1/F2 = 0.51
- Tertia slot F3/F2 = 1.26
- Quinta slot F4/F2 = 1.64
- Nominal slot F5/F2 = 2.19

The paper's summary states that octave, minor-third and fifth relationships characteristic of Western bells were not identified in the examined Russian bells; it also emphasizes complexity of pitch perception and differences across bell size groups. This supports replacing the false universal Western ratios in the Orthodox demo, but does not justify generalizing the Sysoi ratios to every Orthodox bell.

The frequency objective is enabled. Amplitude objective weights are set to zero because the published SPL values are dependent on strike and observation time and cannot be transferred directly as universal modal amplitude ratios. Geometry seed family remains `orthodox`. The optimizer currently applies the target ratios to sorted FEM mode slots 1, 3, 4 and 5 relative to slot 2; these are provisional F1–F5 proxies, not a validated physical identification of Russian bell partial families. The console labels this run F1–F5 to avoid misnaming them as Western Hum/Tierce/Quint/Nominal.

Orthodox evolution budget is now 50 individuals × 50 generations, deterministic seed 20000219, matching the Western demo budget for comparability. This does not itself validate the physical model or guarantee convergence.

## Verification performed

- Linux Release CMake build of `bell`: passed.
- Short Western internal-solver optimization: passed; `best_geometry.json` produced.
- Short Orthodox internal-solver optimization with new target identifier and zero amplitude weights: passed.
- Full Orthodox 440 Hz demonstration at 50 × 50: passed in 29.66 s. Best objective error reported 0.020343; best modal ratios reported Hum 0.571, Tertia 1.374, Quinta 1.638, Nominal 2.169. The nonzero residuals show this is not an exact fit; further optimization/model work remains. The final sorted-mode ratios are only an exploratory mapping to the Sysoi spectrum and require independent mode-family validation.
- Both `synth-geometry` and `synth-contact` generated valid WAV files from the full Orthodox 440 Hz optimized geometry in a Linux smoke run.
- A Windows `.cmd` execution, MSYS2 UCRT64 build, and full CTest suite were not run in this environment.

## Known limitations / next gates

1. The Orthodox acoustic reference currently represents one large Rostov bell only. Obtain measured spectra for multiple Russian bells across mass/size and historical periods before defining a generalized profile.
2. Check whether the optimizer's current mode indexing reliably corresponds to F1–F5 in measured Russian bells. The present objective assumes sorted modes [0,2,3,4] relative to mode 1, as in the existing generic optimizer mapping.
3. The axisymmetric FEM model cannot represent all mode doublets, asymmetry, and non-axisymmetric radiation behavior of actual bells.
4. Compare measured and synthesized spectra/envelopes for reference recordings; do not tune amplitude weights from the single paper's strike-dependent SPL row.
5. Run the six Windows launchers under MSYS2 UCRT64 and inspect output folder isolation, WAV headers/durations, and screen-reader-friendly console logs.

## Sources

- Aldoshina & Nicanorov (2000), AES preprint 5117: https://www.orthodoxbells.com/acoustics/NICANOROV-Acoustics.pdf
- Hibbert et al. (2014), “Partial Frequencies and Chladni’s Law in Church Bells,” Open Journal of Acoustics 4, 70–77, DOI: 10.4236/oja.2014.42007. Its 2,752-bell survey is broad church-bell context, not an Orthodox-only tuning law.
- “Physical modelling techniques for the dynamical characterization and sound synthesis of historical bells” (2021), Heritage Science: https://www.nature.com/articles/s40494-021-00620-2


## Documentation and architecture correction (post-Stage-8 review)

A subsequent source review found that Stage 8 did **not** implement generic arbitrary-length target profiles. `TargetData` still uses four-element arrays and `EvolutionaryOptimizer::computeError` still uses fixed sorted-mode slots. The CLI still checks a hard-coded `tuning_model` allow-list, and currently obtains the geometry seed family from the target JSON `tradition` field. Therefore the phrase “profile is configurable through JSON” applies only to the four legacy ratio/weight values; it does not mean arbitrary target count or arbitrary mode mapping.

The intended selector-based, profile-name-independent architecture and v1 migration contract are now documented in `docs/TARGET_PROFILE_ARCHITECTURE.md`. The implementation remains future work. Stage 8's successful Orthodox run is evidence only for the legacy fixed-slot experimental Sysoi target, not for a generic Orthodox profile engine.
