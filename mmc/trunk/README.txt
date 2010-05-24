===============================================================================
=                       Mesh-based Monte Carlo (MMC)                          =
=                          Multi-threaded Edition                             =
===============================================================================

Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
License: unpublished version
Version: 0.2 (Cheesecake)

-------------------------------------------------------------------------------

Table of Content:

I.  Introduction
II. Download and Compile MMC
III.Running Simulations
IV. Interpreting the Output
V.  Reference

------------------------------------------------------------------------------- 

I.  Introduction

Mesh-based Monte Carlo (MMC) is a 3D Monte Carlo (MC) simulation software 
for photon migration in random turbid media. MMC combines the strengths 
from both MC-based photon migration and finite-element (FE) method: on one 
hand, it can handle low-scattering media as in MC, on the other hand, it 
can use unstructual meshes to represent curved boundary and complex domains,
as in FE. MMC implements a precise ray-tracing process to propagate a photon
using a Plucker coordinate formula. Both the media and the fluence can
be represented by piece-wise-linear basis functions, thus, providing 
additional accuracy. This implementation also supports multi-threaded 
parallel computing and can give a nearly propotional acceleration when
running on multi-core processors.

MMC uses FE meshes to represent complex domains. To generate
an accurate FE mesh for arbitrary object had been a difficult task
in the past. Fortunately, Qianqian along with other developers had 
made great progress to develop a simple-to-use-yet-powerful mesh 
generation tool, iso2mesh [1], which made this task dramatically 
easier. One should also download and install iso2mesh when running 
all the examples from MMC.

We will soon develop a massively-parallel version of MMC by porting
this code to CUDA and OpenCL. This is expected to produce a hundreds
or even thousands fold of acceleration in speed as we had observed
with the GPU-accelerated Monte Carlo code (Monte Carlo eXtreme, or 
MCX [2]), developed by the same author.

-------------------------------------------------------------------------------

II. Download and Compile MMC

The code of MMC is currently developed in the source control system
using Subversion (SVN). To check out the SVN source code, you
should use the following command:

 svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/mmc/trunk mmc

then type the password as "anonymous_user". This will allow you to 
anonymously check out the entire source code tree.

To compile the software, you need to install GNU gcc compiler toolchain
on your system. For Debian/Ubuntu based GNU/Linux systems, you can type

 sudo apt-get install build-essential

and for Fedora/Redhat based GNU/Linux systems, you can type

 sudo yum groupinstall "Development Tools"
 
to install the necessary compilers. To compile the binary supporting
OpenMP multi-threaded computing, your gcc version should be at least 4.2.
To compile the binary supporting SSE4 instructions, gcc version should
be at least 4.3.4. For windows users, you should install MinGW
with a later version of gcc [3]. For Mac OS X users, you can install
Xcode 3 and find gcc or llvm-gcc [4] from the installation.

To compile the program, you should first navigate into the mmc/src folder,
and type

  make release

this will compile a single-threaded optimized binary under mmc/src/bin
folder. Other make options include

  make omp  # this compiles an OpenMP multi-threaded binary
  make prof # this makes a binary to produce profiling info for gprof
  make sse  # this uses SSE4 optimized subroutines for vector oprations
  make      # this produces an non-optimized binary with debugging symbols

If you append "-f makefile_log" at the end of any of the above 
make commands, you will creat a binary named mmc_log, which uses a 
Logistic-Lattice RNG instead of the 48bit POSIX RNG.

You should be able to compile the code with Intel C++ compiler,
AMD C compiler or LLVM. If you see any error message, please 
follow the instruction to fix your compiler settings or install 
the missing libraries.

-------------------------------------------------------------------------------

III.Running Simulations

Before you create/run your own MMC simulations, we suggest you
first going through all the subfolders under the mmc/example folder
and check out the formats of the input files and the scripts for
pre- and post-processings.

Because MMC uses FE mesh in the simulation, you should create
a mesh for your problem domain before you running the simulation.
Fortunately, you can do this fairly straightforwardly using a 
matlab/octave mesh generator, iso2mesh [1], developed by the same 
author. In the mmc/matlab folder, we also provide additional 
functions to generate regular grid-shaped mesh.

The simplest example can be found under the "example/onecube" 
folder. Please run "createmesh" first from matlab/octave to 
create all the mesh files, these include

  elem_onecube.dat    -- tetrahedral element file
  facenb_onecube.dat  -- element neighbors of each face
  node_onecube.dat    -- node coordinates
  prop_onecube.dat    -- optical properties of each element type
  velem_onecube.dat   -- volume of each element

The input file of the example is onecube.inp, where we
specify all the simulation parameters. The mesh files are 
linked through the volume string (specifying the name stub).
To run the simulation, you should run run_test.sh bash
script. If you want to run mmc directly from the command
line, you can do so by typing

../../src/bin/mmc -n 20 -f onecube.inp -s onecube 

where -n specifies the total photon number to be simulated,
-f specifies the input file and -s gives the output file name.
To see all the supported options, run mmc without any parameters.

The above command only runs 20 photons and it will complete
instantly. An output onecube.dat will be saved to record the
normalized (unitary) fluence at each node. If you specify
multiple time-windows from the input file, the output will 
contain multiple blocks with each block corresponding to the
time-domain solution at all nodes computed for each time window.

More sophisticated examples can be found under 
example/validation and example/meshtest folder, where you
can find createmesh script and data analysis script after
you running the simulations.


-------------------------------------------------------------------------------

IV. Interpreting the Output

to be added

-------------------------------------------------------------------------------
V.  Reference

[1] http://iso2mesh.sf.net  -- an image-based surface/volumetric mesh generator
[2] http://mcx.sf.net       -- Monte Carlo eXtreme: a GPU-accelerated MC code
[3] http://sourceforge.net/projects/mingw/files/GCC%20Version%204/
[4] http://developer.apple.com/mac/library/releasenotes/DeveloperTools/RN-llvm-gcc/index.html
