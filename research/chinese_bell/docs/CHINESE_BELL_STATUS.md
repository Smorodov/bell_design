# Chinese-bell research status

**State: RESEARCH_BLOCKED — diagnostic only; not a production optimizer objective.**

## Current implementation boundary

The optional research module reuses the tradition-neutral cylindrical/Fourier FEM foundation. It is detached from the Western optimizer and excluded from normal builds (`BUILD_CHINESE_RESEARCH=OFF`). Localized strike loading uses the actual surface normal and surface metric. Modal association uses mass-weighted subspace information rather than frequency rank alone.

## Durable findings

- Chinese two-tone bells are reported in the literature to produce distinct tones from different strike positions; their non-circular cross-section and paired modal structures are central to the mechanism.
- The structural FEM/eigenpair machinery is numerically tested by the general solver contracts.
- Corrected non-axisymmetric strike integration is numerically consistent for the implemented reduced formulation.
- Modal subspace continuation can be stable, but a stable 2-D subspace is not proof of a physically identified historical mode family.
- The canonical reduced model's candidate pair showed dominant Fourier orders m=1 and m=4, not a robust common-m pair; the meridional nodal topology was not stable enough to infer a physical `(m,n)a/(m,n)b` identity.
- Reciprocal A/B strike selectivity has not been demonstrated robustly.
- Earlier large frequency separations were affected by modal-window/association artifacts and are not accepted as evidence of the historical two-tone mechanism.
- Contact-width values are mathematical spatial scenarios, not identified mallet diameter, contact-patch diameter, FWHM or measured pressure width.

## Decision

Do not tune ellipticity, profile shape, contact width or tracker thresholds to manufacture the expected interval. Do not promote a broad minor/major-third interval band into a sufficient mechanical gate. Do not enable Chinese acoustic radiation/damping objectives or production optimization.

Reopen only when independent modal-shape/geometry data or a clearly falsifiable model extension can distinguish whether the limitation is the reduced structural representation, selected modal family, strike representation or missing reference evidence.

For the acoustic target and literature basis see `CHINESE_BELL_ACOUSTIC_CANON.md` and `CHINESE_BELL_LITERATURE.md`. For the consolidated physical/numerical reasoning see `CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`.
