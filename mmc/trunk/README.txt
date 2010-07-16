===============================================================================
=                       Mesh-based Monte Carlo (MMC)                          =
=                          Multi-threaded Edition                             =
===============================================================================

Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
License: GNU General Public License version 3 (GPL v3), see License.txt
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
parallel computing and can give a nearly proportional acceleration when
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

The details of MMC are reported in the following paper:

*Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing \
in Plücker coordinates," Biomed. Opt. Express 1, 165-175 (2010)

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

 sudo apt-get install gcc

and for Fedora/Redhat based GNU/Linux systems, you can type

 su -c 'yum install gcc'
 
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
  make sse  # this uses SSE4 optimized subroutines for vector operations
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

The full command line options of MMC include the following:
<pre>
usage: mmc <param1> <param2> ...
where possible parameters include (the first item in [] is the default value)
 -i 	       (--interactive) interactive mode
 -f config     (--input)       read config from a file
 -n [0|int]    (--photon)      total photon number
 -b [0|1]      (--reflect)     1 do reflection at internal&external boundaries, 0 no reflection
 -e [0.|float] (--minenergy)   minimum energy level to trigger Russian roulette
 -u [1.|float] (--unitinmm)    define the length unit in mm for the mesh
 -U [1|0]      (--normalize)   1 to normailze the fluence to unitary, 0 to save raw fluence
 -d [1|0]      (--savedet)     1 to save photon info at detectors, 0 not to save
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save
 -s sessionid  (--session)     a string to identify this specific simulation (and output files)
 -h            (--help)        print this message
 -l            (--log)         print messages to a log file instead
 -D [0|int]    (--debug)       print debug information (you can use an integer or
  or                           a string by combining the following debugging flags)
 -D [''|MCBWDIOXATRP]          1 M  photon movement info
                               2 C  print ray-polygon testing details
                               4 B  print Bary centric coordinates
                               8 W  print photon weight changes
                              16 D  print distances
                              32 I  entering a triangle
                              64 O  exiting a triangle
                             128 X  hiting an edge
                             256 A  accumulating weights to the mesh
                             512 T  timing information
                            1024 R  debugging reflection
                            2048 P  show progress bar
       add the numbers together to print mulitple items, or one can use a string
example:
       mmc -n 1000000 -f input.inp -s test -D TP -b 0
</pre>

The simplest example can be found under the "example/onecube" 
folder. Please run "createmesh" first from matlab/octave to 
create all the mesh files, which include

  elem_onecube.dat    -- tetrahedral element file
  facenb_onecube.dat  -- element neighbors of each face
  node_onecube.dat    -- node coordinates
  prop_onecube.dat    -- optical properties of each element type
  velem_onecube.dat   -- volume of each element

The input file of the example is onecube.inp, where we
specify most of the simulation parameters. The input file reads

 100                  # total photon number
 17182818             # RNG seed, negative to generate
 2 8 0.0              # source position (mm)
 0. 0. 1.             # initial incident vector
 0.e+00 5.e-09 5e-10  # time-gates(s): start, end, step
 onecube              # mesh id: name stub to all mesh files
 1 3 10 50            # no longer used
 1 60 10 50           #
 1 60 1  20           #
 1                    #  num of media (not used)
 1.010101 0.01 0.005 1.0  # scat(1/mm), g, mua (1/mm), n
 4       1            # detector number and radius (mm) (not used)
 30.0    20.0    1.0  # detector 1 position (mm)
 30.0    40.0    1.0  # ...
 20.0    30.0    1.0
 40.0    30.0    1.0

The mesh files are linked through the mesh id (specifying the name stub).
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
