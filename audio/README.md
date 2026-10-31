# Reference audio

`bell_round_q100.wav` is the retained circular reference rendering (`q=1`). It is a listening example, not a calibrated sound-pressure recording or an experimental reference measurement.

The adjacent `bell_round_q100._morphing.txt` records how the geometry-listening command generated this file. At `q=1` the coordinate map is the identity. For `q<1`, `synth-geometry` currently uses the documented reduced-order Rayleigh morphing estimate; it must not be interpreted as a fresh direct FEM solution.

To regenerate the standard modal listening sample from `bells/best_bell.dat`, run `make_audio.cmd`.
