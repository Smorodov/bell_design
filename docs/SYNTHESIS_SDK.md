# Bell Synthesis SDK

## Purpose

`BellSynthesis::bell_synthesis` is a standalone C++17 library target for applications that consume precomputed bell modal profiles and synthesize impacts/audio. It is intentionally separate from the optimizer and does not link the FEM solver, mesh generator, optimization objective, CalculiX adapter, or the `bell` command-line executable.

The SDK target currently contains:

- `.bmodal` structural boundary-profile and `.acfield` radiation-cache loading/validation;
- reconstruction of boundary-only modal data;
- internal-clapper geometry and conservative clapper/modal contact integration;
- directional and finite-distance radiation helpers and listener geometry;
- modal pressure evaluation and WAV writing.

The SDK does not calculate a new FEM solution or optimize bell geometry. A profile must be prepared by BellOptimizer or another compatible modal-analysis producer. The host application supplies the expected geometry, material and solver fingerprints when loading the profile.

## CMake consumer

Install BellOptimizer with `cmake --install <build-dir> --prefix <prefix>`. A consumer project can then use:

```cmake
find_package(BellSynthesis CONFIG REQUIRED)
target_link_libraries(my_game PRIVATE BellSynthesis::bell_synthesis)
```

The installed package exports the target, public headers and CMake package version. It requires C++17 and has no third-party runtime dependency. The SDK builds as a static library by default. Set `-DBUILD_SHARED_LIBS=ON` to produce a shared library; the exported CMake target is the same in either case. A stable C ABI is a future integration option, not part of this release.

A minimal consumer is provided in `examples/synthesis_consumer/`. Both static and shared variants were built and linked from an external consumer project on the reference Linux toolchain.

## Runtime separation

The intended application flow is:

1. Offline preparation: optimizer/FEM produces a modal profile bundle.
2. Runtime load: game or host application loads `.bmodal` and `.acfield` and validates fingerprints.
3. Per-event excitation: host supplies clapper parameters and simulates a strike using cached boundary modes.
4. Rendering: host evaluates the precomputed radiation fields for listener/ear geometry and renders samples or writes a WAV.

Strike force, clapper state, contact history, release amplitudes/velocities, listener pose, T60 and output format remain runtime inputs/state; they are not embedded in the reusable profile.

## Current API boundary and limitations

This release establishes a linkable, installable SDK target and public building blocks. It does not yet provide a single high-level `renderStrike()` facade, real-time callback API, streaming audio sink, C ABI, C#/Unity binding, or thread-safe shared profile manager. The sample consumer verifies package discovery, header installation, static linking and profile-loader linkage; it is not a complete game integration demo.

The modal profile remains a boundary-only synthesis profile, not a full FEM restart. The present clapper/contact and exterior-radiation approximations retain their documented physical limitations.
