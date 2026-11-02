# Chinese-bell research handoff — 1.2.126

## Purpose

This is the handoff point for anyone who wants to continue the Chinese two-tone investigation.

The module is intentionally detached from the normal BellOptimizer build. The negative result is part of the deliverable: future work should begin from the documented limitation instead of repeating the same parameter-search path.

## Starting point

The current reduced model successfully provides:

- a validated structural FEM core;
- non-axisymmetric Fourier structural representation;
- corrected localized surface-normal strike integration;
- exact local surface metric for the elliptical mapping;
- modal subspace continuation;
- modal field/topology diagnostics.

The current reduced model does **not** successfully provide a physically identified Chinese paired family.

## The central negative result

The canonical candidate tracked through ellipticity remains numerically continuous as a two-dimensional subspace, but the two members do not exhibit a stable common azimuthal topology:

- subspace similarity: `0.928164`;
- frequencies: approximately `1238.26 Hz` and `1254.26 Hz`;
- dominant Fourier order: `m=1` and `m=4`;
- dominant-harmonic purity: `0.400471` and `0.430326`;
- radial surface zero-crossing counts are probe-dependent rather than a stable common topology.

Therefore the pair cannot be labelled as a physical `(m,n)a/(m,n)b` family from the present reduced FEM.

This is a **model-capability limitation**, not evidence that the structural solver is numerically broken.

## Errors that future work must not repeat

### 1. Do not use frequency rank as family identity

Near-degenerate modes exchange sorted order under geometry changes. A pair selected only by rank can produce an apparently attractive interval while following the wrong physical family.

### 2. Do not accept a large interval as proof of the Chinese mechanism

The earlier 300–370 cent observations were not sufficient evidence. Wider modal windows and continuation exposed branch/window-selection effects.

### 3. Do not infer `(m,n)` from one scalar diagnostic

A dominant Fourier coefficient or one nodal count is insufficient when the surface field is strongly mixed. The topology must be demonstrated across the relevant surface region.

### 4. Do not tune the tracker to force the desired pair

Tracker thresholds must describe evidence, not manufacture it.

### 5. Do not reinterpret mathematical contact widths as physical striker dimensions

The 20/35/50 mm studies are spatial-width scenarios only. No physical mallet/contact law was established.

### 6. Do not compensate unresolved mechanics with acoustic proxies

Radiation, decay, roughness, and loudness cannot establish the missing modal-family identity.

### 7. Do not claim historical authenticity from the reduced model

The literature contract is stronger than the current reduced-model evidence.

## What a successful continuation would have to demonstrate

A future model should be able to identify, independently of frequency rank:

1. a physically coherent `(m,n)` family;
2. two related members `a/b` with complementary nodal/antinodal structure;
3. reproducible correspondence of that structure to the documented strike regions;
4. reciprocal strike selectivity;
5. stability under mesh refinement and reasonable geometric perturbation.

Only after these mechanical conditions are demonstrated should interval, damping, radiation, or acoustic perception become meaningful production candidates.

## Likely next model question

The present evidence points toward asking whether the reduced meridional/Fourier representation is sufficiently expressive for the documented non-axisymmetric paired-mode topology.

If it is not, the correct next step is a more complete non-axisymmetric/3-D structural model, not another sweep of the existing reduced parameters.

## Reuse in the main BellOptimizer

The research did produce techniques that may be useful elsewhere:

- exact surface-metric calculations;
- surface-normal generalized-force integration;
- mass-weighted modal subspace comparison;
- full-field modal diagnostics;
- conservative ambiguity handling.

These should be reviewed against the circular-bell production path separately. They should not be reintroduced merely because they worked in an experimental Chinese context.

## Stereo mode

Stereo/rotational impact-position behavior remains a separate future production task. The Chinese research module is not the place to implement it now.

The production stereo design should be based on the generic circular/elliptic bell excitation and modal pipeline, without importing the unresolved Chinese-family assumptions.
