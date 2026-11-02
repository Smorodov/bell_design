# BellOptimizer 1.2.140 — internal clapper physical model

## Added
- Added `Excitation::simulateInternalClapperContact`: planar rigid physical pendulum coupled to inner-wall Fourier FEM modes through a spherical unilateral Hertz contact.
- The clapper uses explicit mass, pivot inertia, centre-of-mass distance, head-centre distance, head radius, pivot position and initial angular state.
- The default Hertz coefficient is calculated from steel-head and bell elastic properties using the reduced modulus and a locally-flat-wall approximation.
- Equal-and-opposite contact forces couple clapper torque and bell modal generalized forces; the undamped impact tracks total pendulum, modal and contact-potential energy.
- Added `bell synth-clapper`, which renders the post-separation modal state through the outer-wall directional far-field approximation into mono PCM16 WAV.
- Added `clapper_dynamics_contract`, `internal_clapper_audio_smoke`, and `generate_clapper_audio.cmd`.

## Validation
- Full C++17 Release build passed with GCC 14.2.
- Focused clapper/contact/acoustic regression set: 6/6 passed.
- Broad CTest run: 76/77 tests confirmed passed; `cylindrical_unified_ellipticity_audit` exceeded the 300-second execution window and remains unverified (not a reported assertion failure).
- 71-node audio: 6 seconds, mono PCM16, 48 kHz; 32 modes retained.
- First contact 0.075658 s; contact duration 0.001024 s; peak force 22447.3 N; maximum relative conservative energy drift 4.61e-5.

## Limitations
- The 5 kg default clapper mass, dimensions and initial angular velocity are illustrative for the bundled medium-sized geometry, not measured historical hardware data.
- Default pivot inertia is the point-mass-equivalent estimate `m*l_cm^2` unless overridden.
- Hertz contact assumes a locally flat wall; the concave bell's two principal curvatures and elliptical contact patch are not yet modeled.
- Impact is undamped; modal T60 decay begins after first separation. Repeated rebounds, friction, plasticity, contact damping, absolute SPL, room acoustics and stereo are not included.
- The older `synth-contact` remains an external translating-striker diagnostic. `synth-clapper` is the internal-clapper path.

## Impact numerical sweep
- At 31 nodes, 1/2/4 microsecond integration steps gave approximately 22.440 kN peak force and 1.03–1.04 ms contact duration; peak-force spread was below 0.002%.
- At 15/31/71 nodes with 2 microsecond steps, peak force was 22.775/22.440/22.447 kN and contact duration 1.01/1.03/1.02 ms.
- This supports convergence of first-impact metrics for the bundled case only; modal spectral and acoustic waveform convergence remain open.
