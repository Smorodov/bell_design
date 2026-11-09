# Continuation audit — BellOptimizer 1.2.184

Date: 2026-11-08. Source: extracted `BellOptimizer-1.2.184-round-bell-extended-audit-2026-11-08.zip`.

## Reproducible build

- Linux x86_64, GCC 14.2.0, CMake 3.31.6, C++17, Release.
- `BUILD_SMOKE_TESTS=ON`, `BUILD_CHINESE_RESEARCH=OFF`.
- Fresh configure completed successfully.
- Full build completed successfully at 100%; final log: `continuation_full_build_final.log`.
- `CMakeLists.txt` in this audit package contains appended targets for the isolated `research/round_bell_audit` executables compared with `audit_results/CMakeLists_checkpoint.txt`. These are research/test build-target additions; they do not modify optimizer runtime logic. The core production source files were not patched with the candidate-pool experiment.

## CTest

- `ctest -N`: 99 tests discovered.
- Full suite run on this single configuration: **99/99 PASS**, 0 failed.
- Total CTest time: 72.54 s, `-j2`.
- Formerly slow `cylindrical_ellipticity_convergence`: PASS in 39.63 s in the full suite. An isolated rerun also passed in 38.45 s, exit code 0.
- The test exercises a synthetic elliptic patch (`q=0.85`, harmonics {0,1,2,3}); it is not itself a round-bell validation.
- Full log: `continuation_ctest_all.log`; isolated log: `continuation_ellipticity_rerun.log`.

## Existing round-bell research artifacts

- Dataset CSV has 64 rows; 61 valid geometries are represented in the 41/51/61/71/81 modal convergence CSV (2440 mode-pair rows).
- Full-field match CSV has 61 rows; expanded-pool multi-fidelity CSV has 61 rows.
- The nine-geometry multi-harmonic batch is explicitly INCOMPLETE; three individual geometry cases have completed. Do not generalize those three to all 61.
- Thickness full-field mass comparison at 81 nodes has only a three-geometry subset. The broader wall-thickness sweep is a separate diagnostic; do not call it complete full-field convergence across all 61.
- Candidate-pool expansion remains a research patch only. It is not applied to production `src/core/modal_tracking.cpp`.

## Windows portability

NOT RUN. The current execution host is Linux. A reproducible MSYS2 UCRT64 procedure and helper script are included at `docs/WINDOWS_UCRT64_AUDIT.md` and `test_windows_ucrt64.cmd`.

## Integration decision

**NOT READY FOR PRODUCTION INTEGRATION.** The full Linux test suite and build now pass, but Windows UCRT64 remains unverified; full multi-harmonic validation across the geometry set is incomplete; and the candidate-pool patch needs independent review and regression testing before any production adoption. `western_five_partial` and production sources remain unchanged.
