# Bell-bronze modal damping: literature reference

## Purpose

This note establishes a defensible starting value for the **material-only** loss factor used by the internal-clapper listening model. It does not claim that a named alloy composition has a universal damping coefficient, and it does not provide the total decay of a bell in air.

## Convention used by BellOptimizer

For light damping, the free modal amplitude is represented as

\[
A_n(t)=A_n(0)\exp(-\alpha_n t),\qquad
\alpha_n=\frac{\eta_{loss}\omega_n}{2}=\pi\eta_{loss}f_n.
\]

Here `eta_loss` is the structural loss factor. In the conventional critical-damping-ratio notation, `zeta`, the lightly damped relation is approximately `eta_loss = 2*zeta`, and the amplitude decay rate is `alpha = zeta*omega`. This convention conversion matters: values called “damping” in bell papers are not always the same quantity.

## Published bell data

1. P. J. M. Roozen-Kroon, *Structural optimization of bells*, Eindhoven University of Technology PhD thesis, 1992, DOI [10.6100/IR373470](https://doi.org/10.6100/IR373470). Section 2.6 reviews vacuum/low-radiation measurements by Schad (1973) and Van Heuven (1949), and adopts `zeta_m = 0.55e-4` as a provisional bell-bronze material damping value. Later experimental measurements on two carillon bells yielded a mean `zeta_m = 1.37e-4`; the thesis averages the literature and new measurements and uses `zeta_m = 0.9e-4` for subsequent bell analyses. Converted to the loss-factor convention used here, this is `eta_loss ~= 1.8e-4`.
2. M. H. van Houten, A. J. G. Schoofs and D. H. van Campen, “Damping of Bells using Experimental and Numerical Methods,” Fifth International Congress on Sound and Vibration, 1997. Four foundry bells had mean **material** damping ratios of approximately `1.3e-4` to `2.6e-4`; converted loss factors are approximately `2.6e-4` to `5.2e-4`. The paper states that material and acoustic damping are separate contributions and both are typically of order `1e-4`. It also reports that acoustic damping is particularly high for some modes in roughly the 600–1500 Hz range, and that the lowest hum mode can persist longest.

## Current reference preset

`bell_bronze_reference` is assigned `eta_loss=1.8e-4`, the converted 1992 literature mean. This is a **literature reference for bell bronze**, not a composition-specific CuSn grade certification and not a batch-specific casting measurement. The `bronze_aluminum_reference` preset has no damping value assigned because no bell-specific measured value is currently available for it.

The previous exploratory value `eta_loss=0.002` is about 11.1 times the literature reference and about 3.8 times the upper end of the converted four-bell material-only range above. It is therefore unsuitable as a default intrinsic bell-bronze material value, although it can remain useful as an intentionally high-loss sensitivity case.

## Critical limitation: total bell decay

A real bell mode loses energy through at least:

- internal material loss in the bronze;
- acoustic radiation into the surrounding air;
- potentially support/suspension and other boundary losses.

Published bell studies estimate acoustic damping separately and show that it depends strongly on modal shape, frequency and bell geometry. BellOptimizer 1.2.149 still applies the material loss envelope after clapper separation and calculates a directional far-field pressure proxy, but it does **not** yet calculate total radiated acoustic power and convert that power into a modal radiation-damping coefficient. Consequently, a `bell_bronze_reference` preview using structural loss alone is a material-only diagnostic and may ring longer than a real bell in air. This is intentional: it prevents acoustic losses from being hidden inside an alloy parameter.

## Next validation stage

1. Keep geometry, clapper, contact integration, mesh and listener direction fixed.
2. Render a 10-second `eta_loss=1.8e-4` material-only reference and a 10-second `eta_loss=0.002` high-loss comparison.
3. Compare per-mode envelopes and identify which audible partials dominate the perceived tail.
4. Next, estimate per-mode acoustic radiation damping from radiated power over the full far-field sphere and modal mechanical energy. Do not infer it from a single listener-direction pressure value.
5. Validate total modal decay against measured recordings or experimental modal data before changing optimizer objectives.

## Sources

- Roozen-Kroon (1992), [TU/e repository PDF](https://pure.tue.nl/ws/files/2474909/373470.pdf), especially §2.6 and §3.7.
- van Houten, Schoofs & van Campen (1997), [“Damping of Bells using Experimental and Numerical Methods”](https://www.acoustics.asn.au/conference_proceedings/ICSVS-1997/pdf/scan/sv970230.pdf).
