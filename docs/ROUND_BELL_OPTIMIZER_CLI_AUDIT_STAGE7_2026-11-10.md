> Historical Stage 7 note. Its statement that the Orthodox demonstration uses a Western five-partial target was superseded by Stage 8. Current Stage 8 target behavior and limitations are documented in `ROUND_BELL_AUDIO_AND_ORTHODOX_PROFILE_STAGE8_2026-11-10.md` and `TARGET_PROFILE_ARCHITECTURE.md`.

# Round-bell optimizer CLI and launch-script audit — Stage 7

Date: 2026-11-10  
Base: BellOptimizer 1.2.184 Stage 6  
Scope: command-line execution, user-guide correctness, repeatable per-profile/per-frequency launch examples. No production objective or version bump.

## Findings in Stage 6

1. **Relative `--workdir` was applied twice for `optimize`.** `main.cpp` changed the process directory before dispatch, then `command_optimize()` parsed `--workdir` and changed it again. For a relative value such as `bells`, the second change attempted `bells/bells`. This made the documented command unreliable.
2. **The user guide incorrectly said `build.cmd` runs CTest.** The script configures and builds; `test.cmd` is a separate CTest invocation.
3. **The guide omitted the supported `--evolution` and `--targets` options** from the optimize synopsis and did not explain that each independent run needs its own working directory.
4. **Configuration read/parse failures silently continued with defaults.** This could launch an optimization with unintended target/evolution settings after a typo or missing file.
5. **Orthodox naming could be misunderstood.** The Orthodox seed family and geometry assessment are present, but there is no independently validated Orthodox acoustic tuning law. The supplied Orthodox target is a Western five-partial reference template.

## Changes

- `main.cpp` now leaves `--workdir` parsing to `command_optimize()` for the `optimize` subcommand. Other commands retain global working-directory handling.
- Optimize help now advertises `--workdir`, `--evolution`, and `--targets`.
- Missing or invalid target/evolution files now return exit code 2 instead of silently continuing with defaults.
- `prime_freq_target` and `freq_ratios` are checked for finite positive values.
- Corrected the build/test distinction in `docs/USER_GUIDE.md`.
- Added six ready-to-run target presets: Western and Orthodox geometry families at 220, 440, and 880 Hz.
- Added a common Windows launcher and six shortcut `.cmd` files. Every case has an isolated `runs/<profile>_<frequency>Hz` working directory and its own copies of the exact JSON inputs.
- Added `run_all_demo_cases.cmd` for sequential execution and a Russian quick-start guide.

## Verification performed

- Clean Linux CMake Release configuration with `BUILD_SMOKE_TESTS=OFF`.
- `bell` executable compiled successfully (100%).
- `bell --help` displays the updated optimize syntax.
- Relative `--workdir` dispatch check: from `/tmp`, `--workdir bell_workdir_once` reached the requested directory and reported the intentionally missing target file there; it did not attempt `bell_workdir_once/bell_workdir_once`. Expected exit code 2.
- All six target JSON files parsed and passed frequency/ratio/tradition sanity checks.
- Short deterministic Western smoke optimization completed successfully and wrote `best_geometry.json`.
- Short Orthodox smoke optimization completed successfully, selected `Model family: orthodox`, accepted `western_five_partial_reference`, printed the Orthodox geometry assessment, and wrote `best_geometry.json`.

The full Windows MSYS2 UCRT64 build and the complete CTest suite were not run in this environment. The `.cmd` launchers still need the user's Windows-side execution check.

## Interpretation and limits

The Orthodox runs are geometry/style-prior experiments using a Western five-partial reference objective. They are not an Orthodox-specific acoustic-tuning implementation. The Orthodox demo search is 8 individuals × 6 generations; the Western demo search is 50 × 50, so their optimization effort is not directly comparable.

This stage does not change `western_five_partial`, does not enable FEM candidate screening, and does not promote the round-bell modal-family adapter to production. Project version remains 1.2.184.
