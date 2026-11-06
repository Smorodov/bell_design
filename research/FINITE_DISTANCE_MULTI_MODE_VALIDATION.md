# Finite-distance acoustic kernel: multi-mode validation protocol

## Purpose

Establish whether the finite-distance surface integral is numerically stable
for actual cylindrical-FEM modal surface fields over ordinary listening ranges
before considering integration into stereo synthesis.

## Test matrix

- Circumferential content requested from FEM: harmonic orders 0 and 1.
- Receiver distances: 1, 2, 3, and 5 m.
- Both receiver points are evaluated for every returned mode.
- Azimuthal quadrature: 360, 720, and 1440 samples.
- Meridional integration: composite midpoint subdivisions along each existing
  piecewise-linear profile segment; real-FEM gate compares 8 and 16.
- Primary metric: relative complex-field difference between 720 and 1440
  samples, normalized by the magnitude of the 1440-sample result.
- Secondary guard: relative difference between 360 and 720 samples.

Acceptance thresholds in the contract are 0.5% for 720 -> 1440 and 2% for
360 -> 720. These are numerical quadrature gates, not physical-model accuracy
claims. Near a modal radiation null, relative error becomes ill-conditioned;
the test also rejects effectively zero reference fields rather than treating
such a case as a meaningful relative-error result.

## Physical scope

The implemented kernel approximates radiation from the sampled outer shell
using a free-space Green kernel and the normal component of modal surface
displacement. It is not a complete exterior acoustic solution for an unbaffled
open bell. In particular, it omits inner-wall and rim contributions, acoustic
boundary coupling, diffraction by the bell opening, and scattering by the
listener. It uses composite midpoint quadrature on the meridional profile segments. Refinement
linearly interpolates the existing profile and modal field and therefore does not
replace or refine the underlying FEM mesh.

Passing this test supports only the statement that the implemented numerical
integral is stable under the tested azimuthal and meridional quadrature refinement
for the returned FEM modes and receiver geometry. The 1.2.157 real-FEM gate
requires the 8->16 meridional refinement difference to remain below 0.5%; the
observed worst case was 0.241344%. It does not certify the physical
accuracy of the near field or justify replacing the current production renderer.
