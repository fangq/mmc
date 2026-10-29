= Photon reflection at the boundaries of a simple mesh =

== Introduction ==

This example is the same toy problem as examples/onecube - a 10 mm cube split
into 6 tetrahedra, with 20 photons launched from [2,8,0] along [0 0 1] - but
with the boundary reflection turned on (-b 1, using the refractive index
mismatch at the cube surface) and four detectors. It prints the trajectories
of the photons, with their reflections at the boundary, and the exit and
accumulation sites, to check the reflection code visually.

== Steps ==

1. The mesh files for this example have already been generated. If you want
to regenerate them, open MATLAB or Octave and run "createmesh"

2. Run the simulation by typing
  ./run_test.sh

3. Start MATLAB or Octave, and run "plotmmcdebug" to see the photon trajectories
