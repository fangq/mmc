===============================================================================
=                       Mesh-based Monte Carlo (MMC)                          =
=                     Multi-threaded Edition with SSE4                        =
===============================================================================

Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
License: GNU General Public License version 3 (GPL v3), see License.txt
Version: 0.8.pre (Snow cone)

-------------------------------------------------------------------------------

Table of Content:

I.  Introduction
II. Download and Compile MMC
III.Running Simulations
IV. Plotting the Results
V.  Known issues and TODOs
VI. Reference

------------------------------------------------------------------------------- 

I.  Introduction

Mesh-based Monte Carlo (MMC) is a 3D Monte Carlo (MC) simulation software 
for photon transport in complex turbid media. MMC combines the strength
of both MC-based photon migration and finite-element (FE) method: on one 
hand, it can handle low-scattering media as in MC, on the other hand, it 
can use nonstructural meshes to represent curved boundaries and complex 
domains as in FE. MMC implements a precise ray-tracing technique to propagate 
a photon using a fast Plucker-coordinate-based ray-triangle intersection 
test.  Both the media and the fluence can be represented by piece-wise-linear 
basis functions, thus, providing additional accuracy. This implementation 
also supports multi-threaded parallel computing and can give a nearly 
proportional acceleration when running on multi-core processors.

MMC uses FE meshes to represent a complex domain. To generate
an accurate FE mesh for arbitrary object had been a difficult task
in the past. Fortunately, this had been greatly simplified
with the development of a simple-to-use-yet-powerful mesh 
generation tool, iso2mesh [1]. One should download and 
install the latest iso2mesh toolbox when running all the 
build-in examples in MMC.

We will soon develop a massively-parallel version of MMC by porting
this code to CUDA and OpenCL. This is expected to produce a hundreds
or even thousands fold of acceleration in speed as we had observed
with the GPU-accelerated Monte Carlo code (Monte Carlo eXtreme, or 
MCX [2]), developed by the same author.

The details of MMC can be found in the following paper:

  Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing 
  in Plücker coordinates," Biomed. Opt. Express 1, 165-175 (2010)

The author of this paper is greatly appreciated if you cite the 
above paper as reference if you use MMC and related software
in your publication.

-------------------------------------------------------------------------------

II. Download and Compile MMC

The latest release of MMC can be downloaded from the following URL:

  http://mcx.sourceforge.net/cgi-bin/index.cgi?Download

The development branch (not fully tested) of the code can be accessed 
using Subversion (SVN), however this is not encouraged. To 
check out the SVN source code, you should use the following command:

  svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/mmc/trunk mmc

then type the password as "anonymous_user". This will allow you to 
anonymously check out the entire source code tree.

To compile the software, you need to install GNU gcc compiler toolchain
on your system. For Debian/Ubuntu based GNU/Linux systems, you can type

  sudo apt-get install gcc

and for Fedora/Redhat based GNU/Linux systems, you can type

  su -c 'yum install gcc'
 
to install the necessary compilers. To compile the binary supporting
OpenMP multi-threaded computing, your gcc version should be at least 4.0.
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
make commands, you will create an executable named mmc_log, which uses a 
Logistic-Lattice RNG instead of the 48bit POSIX RNG.

You should be able to compile the code with Intel C++ compiler,
AMD C compiler or LLVM. If you see any error message, please 
follow the instruction to fix your compiler settings or install 
the missing libraries.

After compilation, you can add the path to the "mmc" binary (typically
mmc/src/bin) to the search path, so you don't have to type the fully 
path to run it. To do so, you should modify your PATH environment 
variable. Detailed instructions can be found at [5].

-------------------------------------------------------------------------------

III.Running Simulations

Before you create/run your own MMC simulations, we suggest you
first going through all the subfolders under the mmc/example 
directory and check out the formats of the input files and the 
scripts for pre- and post-processing.

Because MMC uses FE mesh in the simulation, you should create
a mesh for your problem domain before you running the simulation.
Fortunately, you can do this fairly straightforwardly using a 
matlab/octave mesh generator, iso2mesh [1], developed by the same 
author. In the mmc/matlab folder, we also provide additional 
functions to generate regular grid-shaped tetrahedral mesh.

It is HIGHLY recommended to use the "savemmcmesh" function under
mmc/matlab folder to save the mesh produced by iso2mesh, as it
performs a number of tests to ensure the consistency of element 
orientations. If you choose not to use savemmcmesh, you 
MUST call "meshreorient" function in iso2mesh for elem/face
to make sure all elements are oriented in the same direction. 
Otherwise, MMC will give incorrect results.

The full command line options of MMC include the following:
<pre>
usage: mmc <param1> <param2> ...
where possible parameters include (the first item in [] is the default value)
 -i 	       (--interactive) interactive mode
 -s sessionid  (--session)     a string to label all output file names
 -f config     (--input)       read config from a file
 -n [0.|float] (--photon)      total photon number
 -b [0|1]      (--reflect)     1 do reflection at int&ext boundaries, 0 no ref.
 -e [0.|float] (--minenergy)   minimum energy level to trigger Russian roulette
 -U [1|0]      (--normalize)   1 to normalize the fluence to unitary,0 save raw
 -d [1|0]      (--savedet)     1 to save photon info at detectors,0 not to save
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save
 -C [1|0]      (--basisorder)  1 piece-wise-linear basis for fluence,0 constant
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh
 -u [1.|float] (--unitinmm)    define the length unit in mm for the mesh
 -h            (--help)        print this message
 -l            (--log)         print messages to a log file instead
 -E [0|int]    (--seed)        set random-number-generator seed
 -M [P|PHBS]   (--method)      choose ray-tracing algorithm (only use 1 letter)
                               P - Plucker-coordinate ray-tracing algorithm
			       H - Havel's SSE4 ray-tracing algorithm
			       B - partial Badouel's method
			       S - branch-less Badouel's method with SSE
 -D [0|int]    (--debug)       print debug information (you can use an integer
  or                           or a string by combining the following flags)
 -D [''|MCBWDIOXATRP]          1 M  photon movement info
                               2 C  print ray-polygon testing details
                               4 B  print Bary-centric coordinates
                               8 W  print photon weight changes
                              16 D  print distances
                              32 I  entering a triangle
                              64 O  exiting a triangle
                             128 X  hitting an edge
                             256 A  accumulating weights to the mesh
                             512 T  timing information
                            1024 R  debugging reflection
                            2048 P  show progress bar
                            4096 E  exit photon info
      add the numbers together to print mulitple items, or one can use a string
example:
       mmc -n 1000000 -f input.inp -s test -b 0 -D TP
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
specify most of the simulation parameters. The input file follows
the same format as in MCX (certain fields are no-longer used).
It looks like the following

 100                  # total photon number (can be overwriten by -n)
 17182818             # RNG seed, negative to generate
 2. 8. 0.             # source position (mm)
 0. 0. 1.             # initial incident vector
 0.e+00 5.e-09 5e-10  # time-gates(s): start, end, step
 onecube              # mesh id: name stub to all mesh files
 3                    # index of element (starting from 1) which encloses the source
 4       1            # detector number and radius (mm) (not used)
 30.0    20.0    1.0  # detector 1 position (mm)
 30.0    40.0    1.0  # ...
 20.0    30.0    1.0
 40.0    30.0    1.0

The mesh files are linked through the mesh id (a name stub) with a 
format of {node|elem|facenb|velem}_meshid.dat. All files must exist.
If the index to the element that enclosing the source is not known,
please use the "tsearchn" function in matlab/octave to find out.
Examples are provided in mmc/examples/meshtest/createmesh.m.

To run the simulation, you should run run_test.sh bash
script. If you want to run mmc directly from the command
line, you can do so by typing

 ../../src/bin/mmc -n 20 -f onecube.inp -s onecube 

where -n specifies the total photon number to be simulated,
-f specifies the input file and -s gives the output file name.
To see all the supported options, run mmc without any parameters.

The above command only runs 20 photons and it will complete
instantly. An output onecube.dat will be saved to record the
normalized (unitary) fluence at each node. If one specifies
multiple time-windows from the input file, the output will 
contain multiple blocks with each block corresponding to the
time-domain solution at all nodes computed for each time window.

More sophisticated examples can be found under 
example/validation and example/meshtest folder, where you
can find createmesh scripts and data analysis script after
you running the simulations.


-------------------------------------------------------------------------------

IV. Plotting the Results

As described above, MMC produces a single output file, named as
"session-id".dat. By default, this file contains the normalized,
i.e. under unitary source, fluence at each node of the mesh. If
multiple time-windows are defined, the output file will contain
multiple blocks of data, with each block being the fluence distribution
at each node at the center point of each time-window. The total
number of blocks equals to the total time-gate number.

To read in the mesh files (tetrahedral elements and nodes), one
can use readmmcnode and readmmcelem function under mmc/matlab
directory. Plotting non-structural meshes in matlab is possible with
interpolation functions such as griddata3. However, it is very
slow for large meshes. In iso2mesh toolbox, a fast mesh slicing
& plotting function, qmeshcut, is very efficient in making 3D
plots of mesh or cross-sections. More details can be found at 
this webpage [6], or "help qmeshcut" in matlab. Another useful
function is plotmesh in iso2mesh toolbox. It has very flexible
syntax to allow users to plot surfaces, volumetric meshes and
cross-section plots. One can use something like

  plotmesh(node,elem,'x<30 & y>30');

to plot a sliced mesh.

Please edit or browse the *.m files under all example subfolder
to find more options to make plot from MMC output.

-------------------------------------------------------------------------------

V. Known issues and TODOs

* MMC only supports linear tetrahedral elements at this point. Quadratic \
 elements will be added later
* currently, this code only support element-based optical properties; \
 nodal-based optical properties (for continuous varying media) will be \
 added in the next release
* the current version of MMC does not support saving partial-path-length \
 data at detector sites as MCX does; this is expected to be added in the \
 next release.

-------------------------------------------------------------------------------
VI.  Reference

[1] http://iso2mesh.sf.net  -- an image-based surface/volumetric mesh generator
[2] http://mcx.sf.net       -- Monte Carlo eXtreme: a GPU-accelerated MC code
[3] http://sourceforge.net/projects/mingw/files/GCC%20Version%204/
[4] http://developer.apple.com/mac/library/releasenotes/DeveloperTools/RN-llvm-gcc/index.html
[5] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/AddPath
[6] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?fun/qmeshcut
