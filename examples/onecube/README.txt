= Visual debugging of photon migration in a simple mesh =

== Introduction ==

In this example, we run photon migration in a toy-problem to show
the basic steps and behaviors of the simulation code.
The mesh used in this case is a simple cube, splitted into
6 tetrahedra. The edge length of the cube is 10mm. 
A total of 20 photons are simulated starting
from a source position [2,8,0] along an incident
vector [0 0 1]. We run the simulation and print out 
the moving trajectories of the photon, and the exit and 
accumulation sites. Using a matlab script, one can visualize
the photon moving from the output of the simulation.

== Steps ==

1. The mesh files for this example have already been generated.
If you want to regenerate the mesh files, please open MATLAB
or Octave and run "createmesh"

2. Run the simulation by typing
  ./run_test.sh

3. Start MATLAB or Octave, run "plotmmcdebug" to see the photon trajectory plot

The script runs on the CPU (-G -1) and prints the moves of the photons (-D M).
The same input is also given as JSON files, which read the same mesh files:
onecube.json, and onecube_isotropic.json for an isotropic source, e.g.

  ../../bin/mmc -f onecube.json -n 20 -G -1 -D M
