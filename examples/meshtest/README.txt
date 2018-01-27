= Validation of MMCM in a heterogeneous domain - a sphere inside a cube =

== Introduction ==

In this example, we validate MMCM algorithm using a sphere object
in a homogeneous cubic background. This simulation generates the results
shown in Fig. 3 in the paper.

The cubic domain has a dimension of 60x60x60 mm with optical properties
mua=0.002, mus=1, n=1.37 and g=0.01. A sphere is centered at [30 30 30]mm
with a radius of 10mm. The optical properties of the sphere is 
mua=0.05, mus=5, n=1.37 and g=0.9. The analytical solution, approximated
by a sphere inside a infinite slab, can be computed by the
sphdiffusionslab function in mmc/matlab/ directory (this has already
been done, the result for plane y=30 is saved in the
sphdiffsemiinf.mat file).

To validate MMCM, we generate an FE-mesh for the sphere and the cube.
Three meshes are tested: mesh0: a coarse FE mesh with only 10,000 nodes,
mesh1: a uniform dense FE mesh with 60000 nodes, and mesh2: a mesh
with higher density around the sphere surface and near the source.
The case mesh1 and mesh2 correspond to "MMCM Mesh 1" and "MMCM Mesh 2"
in the paper (Table 1), respectively.

== Steps ==

1. First, you need to create the 3 meshes to run this example. You 
need to first download and install iso2mesh version 1.0 or newer from 

  http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download

2. Start matlab, run createmesh to create all mesh files

3. Open initial_elem.txt, and type the integer in each line to replace
the second integer at the 7-th line of sph{1,2,3}.inp 
(this sets the initial element ID, because the mesh generator may generate
a different mesh on your system, you have to update this manually)

4. run the simulation bash script run_test.sh

5. when all 3 simulations are complete, you start matlab again, and 
run plotmmcsph to generate the plots.
