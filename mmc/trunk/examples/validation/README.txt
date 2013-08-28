= Validation of MMCM in a homogeneous cubic domain =

== Introduction ==

In this example, we validate the MMCM algorithm using a homogeneous 
cubic domain. This simulation generates the results
shown in Fig. 2 in the paper.

The cubic domain has a dimension of 60x60x60 mm with optical properties
mua=0.001, mus=1, n=1.0 and g=0.01. The analytical solution
can be computed by the cwdiffusion (for CW solution) and 
tddiffusion (for time-domain solutions) functions in the 
package of Monte Carlo eXtreme (MCX) under the mcx/utils directory.


== Steps ==

1. First, you need to create the 3 meshes to run this example. You
need to first download and install iso2mesh version 1.0 or newer from

  http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download

2. Start matlab, run createmesh to create all mesh files. Matlab will
print a value for variable eid. This value indicates the initial element
ID in the mesh that encloses the source. You need to manually set this 
number as the second integer on the 7-th line in cube.inp input file.

3. Run the simulation bash script run_test.sh, this will run 30000000
photons with the specified mesh.

4. To produce Fig. 2, you need to run MCX simulations.
You have to download MCX package from 

 svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/mcextreme/mcextreme_cuda/trunk/ mcx

using the password "anonymous_user". You can also download the pre-compiled
binary package from http://mcx.sf.net/cgi-bin/index.cgi?Download

5. Assuming you have compiled and installed MCX, you can run MCX simulation by

5.1 go to mmc/examples/mcxsph run createmcxbin from matlab
5.2 open benchbox.sh, edit the thread number (-t) and thread block size (-T)
based on the compute capability of your card. For 8800GT/9800GT, the -T number
can not exceed 128; for 280/285/295, -T can not be more than 256; for 470, 
-T can not be more than 576 (for MCX 0.4.9 or newer, using -A option is 
recommended)
5.3 run benchbox.sh to generate the MCX output box.mc2

6. When all simulations are done, you start matlab again, and 
run plotcuberes to generate the plots.
