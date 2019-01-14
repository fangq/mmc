= Experimental: Fluence variation vs photon numbers =

== Introduction ==

== Steps ==

1. First, you need to create the mesh files to run this example. You
need to first download and install iso2mesh version 1.0 or newer from

  http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download

2. Start matlab, run createmesh to create all mesh files. Matlab will
print a value for variable eid. This value indicates the initial element
ID in the mesh that encloses the source. You need to open vartest.json
with a text editor and replace the number after "InitElem".

3. Run the simulation bash script run_test.sh, this will run the simulation
for 10 repetitions and save fluence maps for 1E4, 1E5, 1E6 and 1E7 photons.

4. When all simulations complete, you start matlab again, and 
run plotstdhist.m to generate a plot between signal standard error
and photon number.


*This test is experimental. The interpretations to the results need
further understandings.
