# BellOptimizer 1.2.184 — Stage 17 checkpoint (2026-11-11)

## Included

- Carries forward the Stage 16 round-bell scale/energy investigation and its diagnostics.
- Consolidates the user/developer reference into `docs/PROJECT_GUIDE.md`, `docs/DEVELOPMENT.md`, `docs/VALIDATION.md` and the Russian optimization examples.
- Removes six overlapping/outdated guide files and generated raw modal-overlap CSV result matrices. Keeps the compact numerical summaries and the source geometry dataset needed to reproduce research probes.
- Removes the obsolete Stage 6 checksum list and all local build products from the archive tree.

## Verified

- Clean GCC 14.2 / C++17 Release build: `bell` and `BellSynthesis::bell_synthesis`.
- Eight focused contracts passed: synthesis SDK link, multi-fidelity pipeline, physical modal selection, optimizer shape/scale, hammer contact, clapper dynamics, clapper scaling and clapper scale/energy audit.
- `clapper_scale_energy_audit` confirms its synthetic similarity contract: initial energy follows `s^4`, energy closure is within tolerance, and the synthetic transfer fraction remains within the stated tolerance.

## Important interpretation

The Stage 16 full FEM scale sweep is documented in `docs/VALIDATION.md`. Under its selected geometrically similar clapper setup, larger bells received more absolute modal energy and the transfer fraction rose slightly (13.31% to 15.12%); this does not support the idea that the large bell received less mechanical energy in that setup. However, the strike occurred at 33.9% height, outside the typical 50–70% band, so the setup is not yet representative of a conventional installation.

The sweep also shows that the fixed modal request truncates the high-frequency band increasingly as bell size grows. That is a plausible source of the large bell's duller/perceptually quieter sound. All contact WAV outputs are peak-normalized to 0.95, so their digital peak/RMS cannot be used to compare absolute transferred energy or calibrated SPL.

## Not claimed as complete

- Absolute SPL calibration and measured clapper/bell validation.
- Modal shape/family purity classification.
- Adaptive absolute-frequency modal expansion for every bell scale.
- Representative strike-trajectory calibration and mesh/mode/time-step convergence.
- Multi-fidelity shortlist recall across the full optimizer search domain.
- Full all-target build/CTest and Windows/MSYS2 validation.
