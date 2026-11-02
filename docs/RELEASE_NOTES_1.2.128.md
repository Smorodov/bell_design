# BellOptimizer 1.2.128

## Acoustic-field / stereo architecture audit

This release is a documentation and architecture-gate release. Production acoustic synthesis behavior is intentionally unchanged.

### Findings

- The current scalar production acoustic coefficient is sufficient for monophonic modal synthesis but is not a spatial acoustic field.
- An exactly axisymmetric circular-bell field is the `m=0` harmonic and is invariant under azimuthal rotation.
- Therefore meaningful stereo rotation requires retained non-axisymmetric excitation/radiation information; a scalar amplitude plus a stereo angle would be only panning, not physical field rotation.
- Existing vector-FEM Fourier modal samples remain structural data and are not reclassified as acoustic radiation data.
- The future production target is a compact far-field directional representation with complex azimuthal coefficients and an explicit reference-frame contract.
- Near-field spatial acoustic fields are deliberately out of scope for the first stereo implementation.

### No behavior change

The current scalar synthesis expression and existing optimizer objective remain unchanged. No stereo API has been added.

### Documentation

Updated:

- `docs/ACOUSTIC_FIELD_AND_STEREO_AUDIT.md`
- `docs/ARCHITECTURE.md`
- `docs/DEVELOPMENT_PLAN.md`
- `README.md`

The next implementation gate is the definition and mathematical validation of the compact directional acoustic-field representation.
