# Chinese-bell model status

## Scope

The Chinese-bell model is an isolated **diagnostic / experimental** branch. It is not connected to the production Western optimizer objective.

The historical exploratory sequence through 1.2.121 is consolidated in `CHINESE_TWO_TONE_CONSOLIDATED_AUDIT.md`. The maintained external references are listed in `CHINESE_BELL_LITERATURE.md`.

## Current implementation

- non-axisymmetric bell geometry is represented by the existing cylindrical/Fourier FEM path;
- localized Chinese strike excitation uses the actual surface normal and surface metric during generalized-force assembly;
- modal continuation uses two-dimensional mass-weighted subspace association rather than frequency rank;
- Chinese objective and modal tracker remain test/diagnostic components;
- contact width is a mathematical scenario parameter, not a calibrated physical striker dimension.

## Current scientific status

### Established numerically

- structural FEM matrix/eigenpair contracts remain valid;
- the corrected surface-integrated strike model is numerically consistent;
- modal continuation can remain stable when the spectral window is sufficient;
- profile, mesh, modal association and strike location materially affect the observed response.

### Not established physically

- a robust historical `(m,n)a/(m,n)b` modal family in the current reduced model;
- reciprocal A/B strike selectivity for such a family;
- a physically calibrated striker/contact footprint;
- a validated Chinese radiation/damping model;
- historical authenticity of the reduced model.

### Explicitly rejected shortcuts

- selecting the physical pair by frequency rank alone;
- treating a large frequency split as proof of the two-tone mechanism;
- identifying 20/35/50 mm with mallet diameter, contact-patch diameter or FWHM;
- inserting generic Hertz/contact parameters without independent justification;
- weakening tracker gates to force a preferred branch;
- promoting the minor/major-third interval into a production gate by itself.

## Next stage

The next scientific task is modal-family identification from spatial/topological invariants, followed by A/B strike complementarity. If the current reduced model cannot reproduce that mechanism, the limitation will be recorded rather than compensated with additional tuning parameters.
