# Stage 5 runtime smoke test

This directory contains the exact small-run configuration used to exercise the optional diagnostic path inside `EvolutionaryOptimizer::evaluateCandidate()`.

Run from the project root:

```sh
./build_stage5/bell optimize --solver=internal --workdir stage5_runtime_test
```

The config uses a deterministic seed, population size 2, one generation, and enables `round_fem_diagnostics.csv`. The test output and generated diagnostic CSV are stored at the project root as `audit_results_stage5_runtime.log` and in this directory, respectively.
