# Windows / MSYS2 UCRT64 audit protocol

This protocol is prepared but has **not** been executed in this Linux environment.
It must be run on the user's Windows machine in the MSYS2 UCRT64 environment.

## Preconditions

- Install/update MSYS2 UCRT64 packages: `mingw-w64-ucrt-x86_64-gcc`, `mingw-w64-ucrt-x86_64-cmake`, `mingw-w64-ucrt-x86_64-ninja` (or use MinGW Makefiles), and `make` if needed.
- Open the **UCRT64** terminal, not MSYS, MINGW64, or CLANG64.
- Extract this archive into a writable path without non-ASCII characters for the first run.
- Confirm `echo $MSYSTEM` prints `UCRT64`.

## Commands (UCRT64 bash)

```bash
cd /c/path/to/BellOptimizer
cmake -S . -B build_ucrt64 -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SMOKE_TESTS=ON \
  -DBUILD_CHINESE_RESEARCH=OFF \
  2>&1 | tee audit_results/windows_ucrt64_configure.log
cmake --build build_ucrt64 --parallel 4 \
  2>&1 | tee audit_results/windows_ucrt64_build.log
ctest --test-dir build_ucrt64 --output-on-failure -j2 \
  2>&1 | tee audit_results/windows_ucrt64_ctest.log
sha256sum CMakeLists.txt src/core/modal_tracking.cpp src/core/bell_optimizer.cpp \
  > audit_results/windows_ucrt64_sha256.txt
```

Alternatively, run `test_windows_ucrt64.cmd` from a UCRT64 shell after ensuring
`cmake`, `g++`, and `ctest` are on PATH. If using `cmd.exe`, the shell environment
must be inherited from the UCRT64 terminal.

## Acceptance criteria

1. Configure and build exit code 0.
2. `ctest` reports `100% tests passed, 0 tests failed out of 99`.
3. No missing DLL/runtime dependency when executing the tests.
4. Compare source SHA-256 values against `audit_results/continuation_source_sha256.txt`.
5. Preserve the complete logs and environment file; do not report a partial run as PASS.

If Ninja is unavailable, use `-G "MinGW Makefiles"` and `--parallel 2`.

## Scope boundary

This is a portability check of checkpoint 1.2.184. Do not apply
`research/round_bell_audit/modal_tracking_candidate_pool_expansion.patch` to the
production source while running it. Do not change `western_five_partial`.
