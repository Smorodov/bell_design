# Cross-mesh modal-field comparison

## Purpose and status

`src/fem/cylindrical_modal_field.*` provides an independent diagnostic comparison of complete cylindrical vector-FEM modal fields. It is compiled into the test-only `bell_experimental` library and is not connected to `ModalStateTracker` or the production optimizer.

The input is a non-owning view of a `BellMeshData`, its ordered nonnegative Fourier harmonics, one full global modal vector, ellipticity and material density. The vector layout is the solver layout: `(node, displacement component [ur, uphi, uz], Fourier basis)`, with each nonzero harmonic represented by cosine then sine coefficients.

## Cross-mass product

At a meridional point `(r,z)` and azimuth `phi`, each field is reconstructed with its own Q4 interpolation and Fourier expansion. The cross product is

\[
G_{ab}=\int_{\Omega_{\rm common}}\rho\,\boldsymbol\phi_a\cdot\boldsymbol\phi_b\,dV.
\]

For the elliptic mapping used by the FEM, the physical volume factor in reference meridional coordinates is

\[
dV=\alpha(\phi)^2 r\,|J_{rz}|\,d\xi\,d\eta\,d\phi.
\]

The integration is performed in both directions (quadrature points of mesh A located in mesh B, then mesh B located in mesh A); the two cross integrals are averaged. Each field is evaluated in its own Q4 element. Point location uses element bounding boxes in a uniform spatial bin index followed by inverse bilinear mapping.

The current quadrature uses three-point Gauss-Legendre integration in each meridional coordinate and periodic trapezoidal integration in azimuth. The angular sample count must be sufficient to resolve the product of the two Fourier bases. The default is 64.

For a single pair, the normalized absolute overlap is

\[
S_{ab}=\frac{|G_{ab}|}{\sqrt{G_{aa}G_{bb}}}.
\]

For two-dimensional families, the implementation forms both within-family Gram matrices and the cross-Gram matrix, orthonormalizes each family in its own mass metric, and returns the smallest singular value of the normalized cross-Gram matrix. This is the minimum principal overlap and is invariant to sign changes and rotations of basis vectors within either subspace.

## Domain and material contract

The current implementation requires:

- equal ellipticity for the two fields;
- equal material density (within a small floating-point tolerance);
- valid Q4 connectivity and positive element Jacobians;
- a common physical meridional domain covering at least the caller-specified fraction (default 0.98) of each field's integrated volume.

A low common-domain fraction returns an invalid pair or zero subspace similarity. This prevents a high normalized overlap on a small accidental intersection from being accepted as a valid continuation.

The common-volume fraction is estimated by the same quadrature-point coverage procedure. It is a numerical coverage diagnostic, not an exact polygonal intersection measure; very thin domain differences below quadrature resolution require a dedicated mesh-intersection method.

## Validation evidence

`cylindrical_modal_field_cross_mesh` checks:

- an identical analytic field on 2x2 and 5x4 meshes;
- nested Fourier bases `{0}` and `{0,2}` with the added harmonic absent;
- a two-dimensional family on different meshes;
- rejection of a clipped domain (measured common-volume fraction about 0.8611);
- actual direct-FEM modes from 4x4 and 8x8 meshes, matched by physical full-field overlap rather than sorted mode number.

The original real-FEM contract test found a best cross-mesh modal overlap of 1.0 for its selected low-order mode pair. Version 1.2.92 adds `cylindrical_modal_field_fem_matrix`, which follows the previously identified nearly-degenerate 11th/12th modal family (zero-based indices 10/11) in a production bell geometry at q=0.85. With 14 requested modes, the measured minimum principal overlaps were:

| Continuation change | Minimum 2-D subspace overlap | Common-domain coverage |
|---|---:|---:|
| 11 to 15 profile nodes, same `{0,2,4}` basis | 0.981226 | 1.000000 |
| 11 nodes, `{0,2,4}` to `{0,2,4,6}` | 0.999712 | 1.000000 |
| 11 to 15 nodes and `{0,2,4}` to `{0,2,4,6}` | 0.980768 | 1.000000 |

The legacy boundary matcher did not produce a shape-tracking assignment for the mesh-refinement case (`used_shape_tracking=false`), while the full-field comparison retained the family. These are targeted results for one geometry, material and q; they are not a general convergence proof or evidence that the eigenfrequencies themselves have converged.

## Limitations and next steps

- The two fields must have equal ellipticity. Continuation across changing q needs a reference-configuration pullback or an explicit physical-domain intersection treatment.
- The integral is a consistent continuum mass inner product evaluated with 3x3 meridional quadrature. It is not bitwise identical to the assembled FEM mass matrix, whose element kernel currently uses 2x2 meridional quadrature.
- The implementation supports one- and two-dimensional subspaces, matching the current Chinese-bell pair use case.
- The utility is not yet wired into the Kalman state tracker. The 1.2.92 FEM matrix is a first targeted validation, not yet a wide parameter sweep. Integration should follow additional real-bell studies over mesh refinement and Fourier enrichment.
- The uniform-bin point locator is dependency-free; large production meshes should be benchmarked before using this utility inside an inner optimization loop.


## Version 1.2.92 continuation matrix

`tests/cylindrical_modal_field_fem_matrix.cpp` performs three direct-FEM comparisons on the same bell geometry (q=0.85), density and material. It fixes the reference family to the 11th/12th modes identified by the Chinese-bell readiness study and searches nearby candidate pairs (zero-based indices 8–13). The three scenarios are mesh-only refinement (11 to 15 profile nodes), Fourier-only enrichment (`{0,2,4}` to `{0,2,4,6}`), and simultaneous mesh/Fourier change. Angular integration uses 32 periodic samples, which resolves the maximum sum of retained harmonic orders (6+6).

The measured similarities are 0.981226, 0.999712 and 0.980768 respectively, with full common-domain coverage in all cases. The old boundary matcher did not produce a shape assignment in the mesh-only case. This demonstrates that the cross-mesh field metric can preserve this selected two-dimensional family across the tested discretization changes. It does not yet prove robustness for other geometries, q values, mode families, or actual Chinese-bell strike selectivity.
