# Internal clapper, contact and excitation model

## Implemented path

`synth-clapper` models a suspended rigid pendulum with a spherical soft-steel tip. The tip contacts the inner bell wall; inner-wall Fourier modal displacement is taken from the actual FEM eigenvectors. Equal-and-opposite generalized contact forces couple the clapper and retained bell modes. At first separation, the modal release displacement/velocity is continued with modal decay and radiated from the outer wall.

The older `synth-contact` command is a separate external translating-striker diagnostic. It is not the internal clapper model.

## Contact law and assumptions

The current elastic contact law is spherical Hertz contact:

\[
F = k_H\,\delta^{3/2},\qquad
k_H=\frac{4}{3}E^*\sqrt{R_h},\qquad
E^*=\left[\frac{1-\nu_h^2}{E_h}+\frac{1-\nu_b^2}{E_b}\right]^{-1}.
\]

This uses a locally flat wall approximation. It does not include the two principal curvatures of the concave bell wall, a non-circular Hertz contact ellipse, contact damping, friction, plasticity, wear or a finite-element contact patch. The impact integration is conservative and undamped; modal T60 decay is applied after separation. The current version models the first impact through first separation, not repeated rebounds.

The rigid pendulum requires actual mass, pivot, centre-of-mass distance, pivot inertia, tip-centre distance, tip radius and initial swing conditions for a calibrated instrument. Defaults are illustrative. A default point-mass inertia estimate is not a substitute for CAD or measurement.

## Evidence and limits

- `clapper_dynamics_contract` checks contact, transfer to modal coordinates and conservative energy accounting in its tested setup.
- `internal_clapper_audio_smoke` checks the geometry-to-WAV path on a coarse mesh.
- On the bundled 5 kg demonstrator case, reducing contact integration step from 4 to 1 μs changed peak force by less than 0.002%; contact duration was about 1.03–1.04 ms and conservative energy drift about 3.15e-5.
- Across 15/31/71 profile nodes at 2 μs, peak force was 22.775/22.440/22.447 kN and first-contact time moved from 76.16 ms (15 nodes) to 75.66 ms (71 nodes). This supports convergence of selected impact-event metrics only.

These results do not establish modal spectral convergence, contact-patch accuracy, physical clapper dimensions or agreement with a measured bell recording. WAV output is peak-normalized listening material, not calibrated SPL.

## Optimizer boundary

The clapper acoustic objective is implemented as a tested primitive, but production candidate evaluation still uses the existing Western objective. The optimizer has not been switched to the clapper model. Any future integration needs explicit approval, calibrated physical parameters, mesh/modal convergence, a defined acoustic reference and quantitative acceptance thresholds.

## Deferred physics

Curvature-corrected contact, contact damping, plasticity, friction, wear, flexible clapper modes, repeated strikes, head diffraction, HRTF, room acoustics and absolute acoustic calibration remain out of scope. No 3-D FEM is planned as part of this implementation stage.
