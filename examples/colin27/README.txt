= Simulation with the Colin27 Human Brain Atlas =

== Introduction ==

In this example, we test MMC with a complex human brain atlas: 
the Colin27 atlas. This example is explained in details in the 
original MMC paper [Fang, BOE 1(1),2010]. Running this example 
allows you to reproduce Fig. 4 in the paper.

This example requires the Colin27 atlas mesh. It can be
downloaded at

http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh

In the original MMC paper, the V1 atlas mesh was used.
In this example, we configure to use the V2 atlas mesh.
We believe the difference caused by the mesh versions 
is negligible.


== Steps ==

1. First, you need to download the Colin27 atlas mesh from

   http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
   
   after download, you should extract file "MMC_Collins_Atlas_Mesh_Version_2L.mat"
   and put it under this folder

2. Start matlab, run createmesh to generate all the mesh files.

3. run the simulation bash script run_test.sh

4. The time-resolved fluence/flux maps will be saved in a 
file named brain.dat. You can then load it to make plots. 
You will need the qmeshcut function from iso2mesh toolbox.
Please refer to examples/meshtest/plotmmcsph.m script for 
more details.
