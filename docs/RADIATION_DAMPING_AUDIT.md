# Modal acoustic-radiation damping audit — 1.2.150

## Decision

Do not add a numeric radiation-damping term to `synth-clapper` yet. The current code can evaluate a complex far-field radiation proxy in one selected direction, but it does not yet define a sufficiently complete acoustic radiation problem from which total radiated power can be trusted. In particular, a listener-direction magnitude is not total radiated power and must never be converted directly into a modal loss factor.

The literature bronze preset remains a **material-only** loss reference (`eta_loss = 1.8e-4`). It is not a total damping calibration. The optimizer remains unchanged.

## What the current code actually computes

`src/commands/synthesize_contact.cpp::directionalRadiation()` evaluates, for one frequency and one observation direction, a complex surface integral of the outer-wall normal modal displacement:

`I(n) = integral_S u_n(x) exp(-i k n.x) dS`.

The implementation uses the outer meridional profile, the Fourier-resolved outer modal boundary displacement, midpoint integration along profile segments, and uniform azimuthal samples. The result is then multiplied by modal acceleration in the audio renderer. It is a directional far-field proxy, not an absolute pressure calibration.

Important limitations:

- only the outer bell wall is included in this radiation integral;
- the open mouth, inner wall, rim diffraction and acoustic coupling between inner and outer surfaces are not represented as a complete boundary-value problem;
- the current renderer evaluates one listener direction, not an integral over the full observation sphere;
- air density, radiation impedance and a fully documented pressure normalization are not carried through to an absolute acoustic power;
- the output WAV is peak-normalized, so its amplitude cannot calibrate radiated watts or SPL;
- no measured directivity/power comparison currently validates this approximation.

The current modal eigenvectors are mass-normalized by the internal FEM path (nominal modal mass near 1 kg); any energy calculation must nevertheless use each mode's recorded `modal_mass` and preserve the coordinate normalization consistently.

## Correct energy bookkeeping for a future implementation

For a lightly damped mode with modal mass `M_n`, angular frequency `omega_n`, and peak generalized displacement `Q_n`, the cycle-averaged stored mechanical energy is

`E_n = 1/2 M_n omega_n^2 Q_n^2`.

If a validated acoustic model supplies the cycle-averaged total radiated power `P_rad,n` for that same modal motion, the amplitude-envelope decay rate is

`alpha_rad,n = P_rad,n / (2 E_n)`.

The equivalent small structural loss factor is

`eta_rad,n = 2 alpha_rad,n / omega_n`.

The material-only amplitude decay currently used by BellOptimizer is

`alpha_mat,n = eta_mat omega_n / 2`.

Only after `P_rad,n` is physically defined and independently checked may the renderer combine rates as

`alpha_total,n = alpha_mat,n + alpha_rad,n + alpha_support,n + ...`.

No support-loss term is currently available; it must not be invented. The derivation assumes linear, weakly damped modal motion and cycle-averaged power. The implementation must state its amplitude/peak/RMS convention explicitly and verify dimensions.

## Why a full-sphere integral alone is not enough

A possible far-field route is to evaluate complex pressure over observation directions and integrate acoustic intensity over a sphere:

`P_rad = r^2 integral_Sphere |p(r,n)|^2 / (2 rho c) dOmega`.

This is valid only when `p` is a physically normalized outgoing far-field pressure for the complete relevant radiating boundary and the far-field regime is established. Integrating the current outer-wall proxy over directions would produce a useful *model-relative diagnostic*, but would not repair the omitted mouth/inner-wall/rim physics or supply independent validation. It must not be labelled calibrated radiation damping.

## Required next implementation gates

1. **Acoustic boundary contract:** document the radiating boundary, open-mouth treatment, exterior/interior coupling assumptions, coordinate origin, normal direction, phase convention, and validity range of the far-field approximation.
2. **Power normalization:** derive and unit-test the conversion from modal surface velocity/displacement to complex far-field pressure and acoustic intensity, including air density and sound speed.
3. **Angular integration:** integrate intensity over the full sphere with a documented quadrature; demonstrate convergence when polar and azimuthal orders are refined.
4. **Modal energy:** use the same eigenvector normalization and `modal_mass` as the FEM mode; verify `P_rad/E` is invariant to arbitrary rescaling of the modal coordinate.
5. **Physical checks:** non-negative power; zero for a zero mode; correct frequency/dimension scaling in analytic limiting cases; stable results under angular and surface-mesh refinement.
6. **Independent comparison:** compare predicted directivity and decay against a measured bell recording or a separately validated acoustic model before promoting the result beyond diagnostic status.

Until gates 1–5 pass, retain material-only decay in the listening renderer and report it honestly. Gate 6 is required before describing the resulting total decay as calibrated for a real bell.

## Stereo readiness and sequencing

Stereo does not require calibrated absolute SPL or a complete room model. It does require a stable *complex directional transfer function* for every retained mode, populated from the same modal response and preserving angular phase.

The existing `AcousticDirectionalField` currently supplies only the harmonic representation and rotation operation; it is not populated by the production clapper renderer. The current `directionalRadiation()` computes a single selected direction directly from the FEM surface data, so it is a useful starting kernel but not yet the reusable per-mode field.

Recommended sequence:

1. retain complex directional radiation coefficients for a declared far-field observation regime;
2. validate them against direct evaluation at sampled angles and verify rotation covariance;
3. construct left/right virtual observation directions from one cached modal response and the same field (no second FEM/contact solve);
4. add propagation distance/time-of-arrival and channel gain consistently;
5. later, optionally add HRTF and room reflections as a separate rendering layer.

The first stereo prototype can be developed in parallel with the radiation-power work once the directional field passes its own validation. It must not wait for absolute radiation-power calibration, and it must not be presented as a full binaural/room simulation.

## Status at 1.2.150

- Material damping reference: available as a literature-based starting preset; not casting-specific calibration.
- Directional far-field proxy: implemented for one listener direction in `synth-clapper` and `synth-contact`.
- Total modal acoustic-radiation damping: **not implemented; blocked on the acoustic boundary and power-normalization contract above**.
- Reusable per-mode directional field in production synthesis: **not connected**.
- Stereo renderer: **not implemented**.
- Optimizer integration: **unchanged and intentionally out of scope**.
