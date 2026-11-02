# Chinese modal-family contract

Status: **RESEARCH / DIAGNOSTIC ONLY**

## Purpose

The Chinese two-tone model must identify a physical modal family before any interval,
acoustic, or production objective is allowed to depend on it.

The literature supports three separate observations:

1. Bell modes are characterized by nodal meridians and nodal circles.
2. Chinese two-tone bells exhibit paired modes; a pair is commonly described as
   `(m,n)a` and `(m,n)b`, with the same modal order but different nodal-line
   placement.
3. Different strike positions selectively excite the two members of the pair.

These points are supported by Rossing's review and by experimental work on modern
Chinese two-tone bells; Jing also states that the lens-shaped cross section is central
to the two-tone mechanism while noting that a full 3-D model is needed for complete
acoustic simulation.

## What the current reduced FEM can measure

The current direct cylindrical FEM already exposes the complete vector eigenfield in
Fourier basis form. Therefore the following diagnostics are legitimate:

- dominant azimuthal/Fourier order `m`;
- continuation of a two-dimensional modal subspace across ellipticity or mesh changes;
- modal orthogonality and generalized-force projection;
- relative excitation by explicit `sui`/`gu` strike-position surrogates.

The new `chinese_modal_family_contract` test formalizes these measurements. It does
**not** create a production gate and it does not assign a historical `(m,n)` label from
frequency rank.

## What remains unresolved

The reduced cylindrical field does not yet provide an independently validated,
geometry-aware identification of the meridional nodal-circle count `n` together with
physical `sui`/`gu` nodal-line locations. Consequently:

- dominant `m` is diagnostic, not historical proof;
- subspace similarity is a numerical continuation measure, not proof of a physical pair;
- strike selectivity is reported, not thresholded;
- no interval target is derived from this test;
- no acoustic objective is enabled by this test.

In particular, an apparently stable two-dimensional eigenspace is insufficient. The
physical contract requires the modal topology and the strike locations to agree with
the documented paired-mode mechanism.

## Decision rule

If a future implementation cannot establish `(m,n)` pair identity and complementary
nodal structure in a way that is stable under mesh refinement and geometry continuation,
the result must be recorded as a model limitation. Tracker tuning, frequency-rank
selection, or contact-width tuning must not be used as substitutes.

## External basis

- M. Jing, *A theoretical study of the vibration and acoustics of ancient Chinese bells*,
  JASA 114(3), 1622–1628 (2003), DOI 10.1121/1.1600727.
- T. D. Rossing, *Acoustics of Eastern and Western bells, old and new*, J. Acoust. Soc.
  Jpn. (E) 10(5), 241–252 (1989).
- J. Tsai, Z. Jiang, T. D. Rossing, *Vibrational Modes of a Modern Chinese Two-Tone Bell*,
  Acoustics Australia 20(1), 17–19 (1992).
- T. Chen, D. Zheng, *Acoustical Properties of Chime Bells*, Chinese Journal of Acoustics
  4(1), 1–10 (1985), DOI 10.15949/j.cnki.0217-9776.1985.01.001.

This document deliberately separates literature-supported physical statements from
measurements that are only established numerically inside BellOptimizer.
