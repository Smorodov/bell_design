# Reproducing the round-bell audit experiments

Run from the project root in a C++17-capable CMake environment.

```sh
cmake -S . -B build_mesh_repeat \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SMOKE_TESTS=ON \
  -DBUILD_CHINESE_RESEARCH=OFF
cmake --build build_mesh_repeat --parallel 4
```

The four principal geometry/modal research executables are:

```sh
./build_mesh_repeat/round_geometry_dataset build_mesh_repeat/round_geometry_dataset.csv
./build_mesh_repeat/round_geometry_modal_convergence \
  build_mesh_repeat/round_geometry_dataset.csv \
  build_mesh_repeat/round_geometry_modal_convergence.csv
./build_mesh_repeat/round_geometry_multifidelity_probe \
  build_mesh_repeat/round_geometry_dataset.csv \
  build_mesh_repeat/round_geometry_multifidelity.csv
./build_mesh_repeat/round_geometry_fullfield_match \
  build_mesh_repeat/round_geometry_dataset.csv \
  build_mesh_repeat/round_geometry_fullfield_match_16.csv 15
```

Thickness sweeps:

```sh
./build_mesh_repeat/round_geometry_thickness_sweep \
  build_mesh_repeat/round_geometry_dataset.csv \
  build_mesh_repeat/round_geometry_thickness_41.csv 41
./build_mesh_repeat/round_geometry_thickness_mass_sweep \
  build_mesh_repeat/round_geometry_dataset_41_subset.csv \
  build_mesh_repeat/round_geometry_thickness_mass_41_3.csv 41
```

The temporary candidate-pool patch is in
`research/round_bell_audit/modal_tracking_candidate_pool_expansion.patch`.
It is **not applied** in the packaged source tree. Apply it only in a separate
working copy when reproducing the expanded-pool probe.

The multi-harmonic probe is exploratory and expensive. It uses 96 direct
candidates and is not part of the round-bell production readiness gate.

The CTest suite contains 99 tests. The recorded completed filtered run passed
98/98 after excluding `cylindrical_ellipticity_convergence`; that one test is
still unresolved and must not be counted as passed. Windows/MSYS2 UCRT64 has
not yet been tested.
