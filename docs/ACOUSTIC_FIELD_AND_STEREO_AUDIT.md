# Acoustic field and stereo — current engineering record

## Status at 1.2.166

The `synth-clapper --stereo` workflow has a shared mechanical/modal history and two receiver evaluations. It is a useful geometric stereo listening prototype, not a calibrated acoustic prediction. The far-field receiver path remains the default; the finite-distance path remains opt-in.

### Shared computation contract

For a given bell, clapper state and strike:

1. solve the structural modes once;
2. simulate pendulum/contact once;
3. derive one post-separation modal state and decay history;
4. evaluate that same history at left and right receiver positions;
5. encode the two channels.

Per-channel FEM, contact integration or modal time evolution is forbidden. Receiver geometry may change without altering the structural/contact result.

### Receiver paths

**Far-field (default).** A complex azimuthal directional field is built per mode and sampled in the two receiver directions. This is the intended reusable representation for head rotation. It remains an approximate free-space radiation model.

**Finite-distance (opt-in).** An outer-shell Green-kernel surface integral is evaluated at the exact ear points. It is a comparative receiver kernel, not a complete exterior acoustic boundary-value solution. Its extra receiver integrals do not repeat FEM/contact/modal dynamics.

Neither path includes head diffraction, pinna/HRTF, torso scattering, room reflections or calibrated sound pressure level. WAVs are peak-normalized listening previews.


## 1.2.166 — shared far-field coefficient construction

The former command-local direct azimuthal integral and Fourier/Bessel field builder have been moved into `src/acoustics/directional_radiation.*`. Production mono (`synth-contact`) and stereo/clapper far-field synthesis now call the same `buildDirectionalRadiationField` implementation. The direct azimuthal quadrature remains available as an independent reference evaluator and is no longer used in production rendering.

A maintained equivalence contract compares complex values (therefore both amplitude and phase) for a synthetic modal surface field over 3 frequencies (320, 1873.25 and 4200 Hz), 3 listener polar angles (0.19, 0.71 and 1.31 rad), and 17 azimuths at each condition. With 8192 direct azimuth samples, the observed maximum relative complex difference was approximately `2.99e-14` (maximum absolute difference approximately `2.02e-18`). This is numerical agreement for the tested field and implementation conventions, not a physical accuracy claim.

The harmonic builder is now the single production far-field path. The direct integral is retained only as a regression oracle so future changes to Fourier sign, phase convention or normalization cannot silently alter the rendered field. The finite-distance Green-kernel path remains separate because it represents a different receiver approximation.

## Consolidated numerical evidence

The following are reproducible engineering comparisons from the archived one-off audits. Raw WAV/CSV outputs have been removed from the source package; the values and scope are retained here and in `research/CONSOLIDATED_NUMERICAL_RESULTS.md`.

### 1.2.161 — real-FEM receiver comparison

- The comparison used a real cylindrical FEM solution, 13 modal components, two ears, listener distances 1/2/3/5 m and head orientations 0°/15°.
- Global absolute complex difference normalized by the strongest modal transfer was approximately 1.18%.
- Relative errors for individual weakly radiating modes could reach about 32%; these are sensitive to near-zero radiation denominators and must not be interpreted as the overall field error.


### 1.2.167 — coordinate covariance and no-copy field reuse

`AcousticDirectionalField::evaluateRotated(theta, frameRotation)` now evaluates a rotated coordinate frame directly from the original Fourier coefficients. It is algebraically equivalent to materializing `field.rotated(frameRotation)` and then evaluating the result, but does not copy or mutate the coefficient vector. This makes repeated receiver-angle queries independent of the number of head-yaw samples once a per-mode field has been built.

The rotation contract now checks: identity/full-turn behavior, the complex phase factor for each azimuthal order, coefficient-norm preservation, additive composition, joint source/listener coordinate covariance, equivalence of no-copy and materialized rotations over multiple angles, and immutability of the base field. These are algebraic tests of the harmonic representation. They do not yet prove that a rotated source, strike excitation, and receiver geometry are wired consistently through every production command.

The far-field stereo renderer already builds one field per mode and evaluates both ears from it; it does not rebuild FEM or the surface integral for the second ear. The new no-copy API extends this to repeated coordinate-frame queries. Exact finite-distance receiver evaluation remains separate. Physical head diffraction remains out of scope; ear-specific polar-angle interpolation is now implemented and has passed the real-FEM numerical contract described in 1.2.171.

### 1.2.162 — complete stereo waveform comparison

Four matched 10-second, 48 kHz stereo previews were rendered from identical FEM/contact/modal release histories (0° and +15° head yaw; both receiver paths).

- At 0° yaw, left/right channels matched to PCM quantization precision.
- At +15° yaw, L/R correlations were 0.9259 (far-field) and 0.9441 (finite-distance).
- Interaural propagation offset was 68.8 μs, approximately 3.30 samples at 48 kHz.
- The finite-distance vs far-field waveform correlations after shared peak normalization were 0.9489 at 0° and 0.9569/0.9402 for left/right at +15°.
- Band-level differences were measurable, but neither path has an independent physical reference that would justify declaring it more accurate.

### 1.2.163 — listener matrix

A real-FEM comparison covered 27 configurations: distance 1/3/5 m, polar angle 0.45/0.82/1.20 rad, yaw 0°/15°/30°, both ears and 11 returned modes (594 receiver/mode samples).

- Aggregate normalized complex L2 difference: **0.920958%**.
- Maximum eligible per-mode/receiver relative difference: **11.49% at 1 m**, **3.78% at 3 m**, **2.33% at 5 m**.
- Individual relative comparisons exclude far-field amplitudes below 0.1% of the strongest mode in the same configuration to avoid near-zero denominators.

These are differences between two approximate kernels, not physical error bounds. The matrix uses one representative bell, fixed 0.18 m ear spacing and idealized receivers.

## Decisions

- Keep far-field as the default; finite-distance stays opt-in.
- Do not claim physical superiority from this comparison alone.
- Do not add HRTF/head diffraction or room acoustics until the free-field geometric layer has passed rotation-covariance and phase/time-convention tests.
- Cache/reuse the immutable per-mode directional field for arbitrary receiver yaw; rotation must not rebuild the FEM solution or radiation field.
- Preserve absolute phase through modal summation. Peak normalization is a listening-output convenience, not SPL calibration.

## Remaining acceptance gates (current baseline 1.2.172)

1. Stereo phase/time convention: verify interaural delay, modal complex phase and sample-rate conversion without double-counting propagation delay.
2. Stereo listening validation: compare generated recordings with controlled measurements only after geometry, material, clapper and microphone placement are known.

The direct-vs-harmonic complex equivalence and coordinate/rotation contracts validate implementation consistency; they do not establish physical accuracy of the radiation approximation.

## 1.2.168 — joint source/strike/listener coordinate contract

Added `stereo_coordinate_covariance_contract`. It rotates the listener centre and both idealized ear points about the bell axis while retaining range, polar angle, ear spacing, yaw relative to the local radial/tangential frame, and propagation delays. It also verifies that the strike-to-listener relative azimuth is unchanged when source, strike and listener coordinates receive the same global rotation. The complex left/right transfer is checked using the no-copy field-frame evaluation.

This closes a stronger *geometry/representation* covariance gate than the isolated Fourier-field test. It is still not a production end-to-end clapper test: the current command does not expose a global bell/source orientation transform, and the test does not rerun contact dynamics under a rotated strike. In particular, FEM modal basis phase conventions and the mapping of `strike_azimuth_rad` into modal release coordinates still need a dedicated integration test before claiming full physical source/strike covariance.

The current far-field stereo renderer samples an azimuth-only field built at the listener-centre polar angle. The two ears have slightly different polar angles at finite range. This is not guaranteed to be a small error: the 1.2.169 synthetic multi-harmonic audit measures the effect explicitly below. A true far-field 3D solution requires polar-angle dependence (for example a spherical-harmonic representation or validated interpolation over polar angle), not merely an azimuthal rotation. The finite-distance receiver path evaluates each ear point independently, but it remains a separate approximate kernel.


## 1.2.169 — clapper/modal covariance and ear polar-angle error

### End-to-end clapper/modal frame covariance

Added `clapper_stereo_frame_covariance_contract`. It constructs an exactly degenerate cosine/sine modal pair, applies an orthogonal rotation to the real Fourier modal basis, rotates `strike_azimuth_rad` by the same angle, and runs the production `simulateInternalClapperContact` integrator for both frames. The test compares first-contact and separation times, peak force, contact point, the complete contact-force history, and modal displacement/velocity at release. It then builds the corresponding production Fourier/Bessel radiation fields and compares the complex modal pressure after rotating the listener azimuth.

This verifies the contact-to-release-to-radiation coordinate contract for an exactly degenerate harmonic pair and the actual clapper dynamics implementation. It does not prove covariance for arbitrary non-degenerate FEM modes, a non-axisymmetric bell geometry, or arbitrary eigensolver basis mixing. In particular, a rotated pair is physically admissible as a basis change only within a degenerate (or numerically unresolved) modal subspace.

### Ear polar-angle approximation measurement

Added `stereo_polar_angle_approximation_contract`. It compares (a) the production shortcut that evaluates both ears using the listener-centre polar angle against (b) direct azimuthal quadrature evaluated at each ear's own polar/azimuth angles. The polar-specific Fourier/Bessel builder is independently checked against direct quadrature at each ear.

For the documented synthetic multi-harmonic fixture (640 Hz, 0.18 m ear spacing, centre polar angle 0.83 rad, fixed yaw and azimuth), the maximum absolute complex error normalized by the strongest direct ear response in the distance sweep was:

| Listener distance | Centre-polar approximation error / reference peak |
|---:|---:|
| 0.75 m | 61.22% |
| 1.5 m | 29.30% |
| 3.0 m | 14.20% |
| 6.0 m | 6.97% |

For this fixture the polar-specific Fourier/Bessel result agreed with direct quadrature to approximately `1.5e-14`–`1.7e-14` of the reference peak. These figures isolate the polar-angle shortcut; they are not a universal bound for real bells and depend on frequency, bell size, modal harmonics, receiver direction and head geometry. The observed error is large enough that the centre-polar shortcut must not be presented as generally accurate for ordinary listening distances.

**Implementation decision:** do not silently change the production renderer in this audit. Keep the current shortcut explicitly labelled as an approximation until a polar-dependent far-field representation is implemented and validated against direct per-ear evaluation on real FEM modes. The next implementation gate is a 3D polar/azimuth field or a controlled polar interpolation with quantified error and a listener-only update path that reuses the modal/contact history.


## 1.2.170 — polar-dependent far-field stereo field

The far-field stereo path no longer samples both ears from one field built at the listener-centre polar angle. For each retained mode the production command builds `AcousticPolarDirectionalField` over the local polar interval covering all supported head-yaw values (-90 to +90 degrees), with no more than 2 degrees between nodes. Each node stores the existing complex azimuthal Fourier/Bessel field. A receiver query linearly interpolates the complex field value between adjacent polar nodes and then evaluates the requested azimuth; equivalently, the complex coefficients can be interpolated because the azimuthal basis/order set is shared across nodes. A global 91-node 0..pi grid remains available for standalone validation. This preserves phase information and makes azimuth a phase-only query, while head yaw changes the ear polar/azimuth query within the precomputed envelope; no surface integral is rebuilt.

The stereo command builds this polar family once per mode and reuses it for both ears. FEM, clapper contact, modal release, and time-domain modal evolution remain shared and are not repeated per channel. Mono keeps the previous single-polar builder. The finite-distance path remains an opt-in direct receiver-point kernel and is unchanged.

The new `polar_directional_field_interpolation_contract` compares interpolated results with independently built exact-polar Fourier/Bessel fields at 8 polar angles and 5 azimuths (40 complex comparisons). For its synthetic multi-harmonic fixture, the maximum interpolation error is 0.442% of the largest reference magnitude. The test also checks that polar interpolation commutes with azimuthal frame rotation. This is a regression contract for the interpolation implementation, not a universal acoustic accuracy bound. Further validation must use real FEM modal fields, several frequencies and geometries, and mesh/profile refinement. The local production grid typically needs only 2–4 surface-integral builds per mode at ordinary listening distances. The synthetic contract checks five yaw values across four receiver distances plus intermediate polar queries; the observed maximum is 0.365% of the fixture reference peak. The global 91-node grid is retained as a reference for interpolation tests, not used by the production renderer.



## 1.2.171 — real-FEM polar interpolation validation

Added `real_fem_polar_interpolation_contract`. Unlike the synthetic 1.2.170 fixture, this test obtains modal surface fields from the production cylindrical FEM modal evaluator and checks the production local polar interpolation against independently rebuilt exact-polar Fourier/Bessel fields. It covers three FEM cases: the baseline profile at 17 nodes, the same profile refined to 25 nodes, and a second bell profile at 17 nodes. Across the cases it exercised 54 returned modes spanning 79.55–735.10 Hz, 1/3/5 m receiver distances, five head-yaw values from -90 to +90 degrees, both ears, and 1,620 complex receiver/mode comparisons. There were also 162 independent direct-azimuthal-quadrature spot checks.

Worst per-mode peak-normalized complex interpolation error was `0.000118175` (0.01182%); maximum phase error was `0.00126688` degrees where the exact reference magnitude was at least 10% of that mode's directional peak. Both are below the predeclared gates of 1% complex error and 1 degree phase error. Exact-polar Fourier/Bessel values agreed with the independent direct quadrature to relative error below `2e-5` in all spot checks. These are numerical representation/interpolation results on the tested FEM profiles and mesh levels; they do not validate the far-field physical radiation approximation, exterior boundary conditions, HRTF, room acoustics, or measured sound pressure. The profile variation is a second representative test geometry, not an exhaustive shape family.

**Decision:** retain the linear complex interpolation and 2-degree maximum polar spacing. No tighter grid or higher-order interpolation is justified by these measured cases. The next gate is a reusable receiver-query API that separates listener orientation from the prepared modal/polar state; this was implemented and validated in 1.2.172.


## 1.2.172 — reusable prepared stereo receiver state

Added `PreparedStereoRadiation`, an immutable far-field receiver-query object. Its preparation step accepts already computed modal surface fields and the bell profile, builds one local polar family per supplied mode over a declared head-yaw envelope, and records the polar-family build count. Repeated `query(headYaw, centreAzimuthOffset)` calls calculate the idealized ear geometry and evaluate left/right complex transfers, range scales and propagation delays. Both head yaw and centre azimuth can change within the prepared polar envelope without rebuilding the field family. They do not have access to the FEM solver, clapper integrator or modal-release calculation; those remain upstream inputs and are not part of the receiver-query operation. Out-of-envelope yaw is rejected rather than silently extrapolated.

The production `synth-clapper --stereo` far-field path now prepares this object after the FEM solve, clapper contact and modal release, then queries it for the requested yaw. The object is also reusable by a longer-lived host that wants to change listener yaw repeatedly without rerunning the upstream simulation. The current CLI still performs one render/query per process invocation; it is not an interactive listener-positioning UI and does not persist a cache between separate invocations. Finite-distance and mono branches remain separate.

`real_fem_polar_interpolation_contract` now also exercises the reusable object for all three real-FEM cases. It makes five yaw values (-90, -45, 0, 45, 90 degrees) at three centre-azimuth offsets (0, +0.37, -1.1 rad), compares both ears with independently rebuilt exact-polar fields, checks finite range/delay values, verifies the prepared mode count remains stable, verifies the polar-field build count is unchanged after every query, and requires rejection outside the declared yaw envelope. The existing interpolation gates remain 1% peak-normalized complex error and 1 degree phase error above 10% modal peak.

This closes the **receiver-query reuse API** gate. It did not by itself provide persistent cross-process caching or measured physical-acoustic accuracy. The phase/time gate was subsequently closed in 1.2.173; see the next section.

## 1.2.173 — reusable geometry-dependent acoustic field cache and phase/time kernel

Added a versioned binary `AcousticFieldCache` artifact. It stores per-mode frequency and the complex Fourier/Bessel directional radiation coefficients on a canonical full polar grid (0..pi). It deliberately excludes excitation-dependent quantities: clapper/hammer force, contact trajectory/history, modal displacement and velocity at release, T60/decay choice, listener ear geometry and PCM samples. A single field cache can therefore be reused for multiple strikes, clapper settings, strike energies, listener distances, centre directions and head-yaw values, provided the bell geometry, material preset and FEM/modal solver configuration are identical.

The cache identity has three required caller-supplied fingerprints: canonical geometry, material preset and FEM/modal solver and radiation-model configuration (including sound-speed/radiation constants). Loading rejects identity mismatches. The file uses an explicit little-endian encoding, a format version, bounded dimensions/strings/payload, finite-value and structural validation, and an FNV-1a corruption checksum. This checksum detects accidental corruption; it is not a cryptographic authenticity mechanism. The `prepareStereoRadiationFromCache` path builds a receiver-query object without rebuilding any radiation surface integral. The cache is currently a library API, not yet wired into CLI flags or an automatic persistent cache manager. Callers must generate stable fingerprints correctly; a wrong fingerprint can defeat cache invalidation.

The canonical full-sphere grid costs more to build and store than the local 2-degree yaw envelope used by one-shot CLI synthesis. This is intentional for a reusable asset that must cover arbitrary receiver positions. The default 91-node grid is a starting policy, not a universal accuracy guarantee; future cache production should carry an interpolation-quality report and permit a denser grid when a geometry/frequency family requires it. Far-field field coefficients are not a BEM/HRTF/room model and the cache does not improve that physical approximation.

Added `stereo_phase_delay_contract` and factored modal pressure evaluation into `Audio::evaluateModalPressureAt`, now used by mono and stereo production synthesis. The signed propagation delay relative to the listener-centre reference is applied exactly once as `t - releaseTime - earDelay` (an ear closer than the centre has a negative relative delay); no additional `exp(-i*omega*delay)` multiplier is applied. The complex directional transfer retains the modal radiation phase, while path length is represented as a time shift. Ear range scaling remains a separate amplitude factor. The contract checks geometric left/right delays, ITD, onset causality, and equivalence to a pure once-only time shift. PCM onset is quantized to the sample grid, so the residual ITD quantization is bounded by one sample.

The cache contract performs a binary round-trip, verifies that the same field can be prepared for three distinct listener configurations and five yaw values without any new surface integrations, and rejects wrong geometry identity and checksum corruption. The test uses a controlled synthetic modal fixture; real-FEM cache-file size/performance and cross-platform compiler interoperability remain follow-up validation items.

## Reusable bell modal profile (1.2.174)

A reusable profile now separates geometry/material/solver-dependent data from strike-dependent state. The logical bundle consists of two versioned, checksummed files:

- `.bmodal`: inner and outer meridional profiles; modal number, eigenvalue, frequency, modal mass and decay-time metadata; inner and outer Fourier boundary traces (radial, azimuthal and axial components plus harmonic labels).
- `.acfield`: per-mode complex directional radiation fields on the canonical polar grid, using the existing acoustic-field cache format.

The bundle is bound to caller-supplied geometry, material-preset and FEM/modal-solver fingerprints. Loading rejects identity mismatches, truncation, checksum errors, inconsistent mode counts and frequency mismatches. The bundle deliberately excludes impact force/energy, clapper trajectory, contact history, modal release coordinates, T60 overrides and PCM output.

This is a boundary-modal synthesis profile, **not** a complete FEM restart: stiffness/mass matrices, interior eigenvectors and optimizer state are not serialized. It is sufficient to reconstruct the boundary `ModalMode` data consumed by the current clapper-contact and directional-radiation paths, provided the caller also supplies the matching profile identity and chooses new strike parameters. A real-FEM round-trip was added to `real_fem_polar_interpolation_contract`: four actual modes from the 17-node baseline case are serialized and restored, and the inner/outer Fourier boundary traces are compared exactly. The restored radiation half can prepare a receiver without rebuilding surface integrals. The CLI has not yet been switched to load this bundle; integration into the fast synthesis command is a separate next step. The two-file bundle is a logical pair rather than a single atomic filesystem transaction; applications should keep both files together and treat a missing half as an invalid profile.
