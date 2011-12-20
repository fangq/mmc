===============================================================================
=                       Mesh-based Monte Carlo (MMC)                          =
=                     Multi-threaded Edition with SSE4                        =
===============================================================================

Author:  Qianqian Fang <fangq at nmr.mgh.harvard.edu>
License: GNU General Public License version 3 (GPL v3), see License.txt
Version: 0.9.0 (Banana Pudding)
URL:     http://mcx.sf.net/mmc/

-------------------------------------------------------------------------------

Table of Content:

I.  Introduction
II. Downloading and Compiling MMC
III.Running Simulations
IV. Plotting the Results
V.  Known Issues and TODOs
VI  Getting Involved
VII.Reference

------------------------------------------------------------------------------- 

I.  Introduction

Mesh-based Monte Carlo (MMC) is a 3D Monte Carlo (MC) simulation software 
for photon transport in complex turbid media. MMC combines the strengths
of the MC-based technique and the finite-element (FE) method: on the 
one hand, it can handle general media, including low-scattering ones, 
as in the MC method; on the other hand, it can use an FE-like tetrahedral 
mesh to represent curved boundaries and complex structures, making it
even more accurate, flexible, and memory efficient. MMC uses the
state-of-the-art ray-tracing techniques to simulate photon propagation in 
a mesh space. It has been extensively optimized for excellent computational
efficiency and portability. MMC currently supports both multi-threaded 
parallel computing and Single Instruction Multiple Data (SIMD) parallism 
to maximize performance on a multi-core processor.

To run an MMC simulation, one has to prepare an FE mesh first to
discretize the problem domain. Image-based 3D mesh generation has been 
a very challenging task only until recently. One can now use a powerful 
yet easy-to-use mesh generator, iso2mesh [1], to make tetrahedral meshes
directly from volumetric medical images. You should download and install 
the latest iso2mesh toolbox in order to run the build-in examples in MMC.

We are working on a massively-parallel version of MMC by porting
this code to CUDA and OpenCL. This is expected to produce a hundred-
or even thousand-fold acceleration in speed similar to what we 
have observed in our GPU-accelerated Monte Carlo software (Monte Carlo 
eXtreme, or MCX [2]).

Please keep in mind that MMC is only a partial implementation of the 
general Mesh-based Monte Carlo Method (MMCM). The limitations and issues
you observed in the current software will likely be removed in the future
version of the software. If you plan to perform comparison studies with 
other works, please communicate with the software author to make
sure you have correctly understood the details of the implementation.
The details of MMCM can be found in the following paper:

  Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing 
  in Plücker coordinates," Biomed. Opt. Express 1, 165-175 (2010)

The author of this paper is greatly appreciated if you can cite 
the above paper as reference if you use MMC and related software
in your publication.

-------------------------------------------------------------------------------

II. Download and Compile MMC

The latest release of MMC can be downloaded from the following URL:

  http://mcx.sourceforge.net/cgi-bin/index.cgi?Download

The development branch (not fully tested) of the code can be accessed 
using Subversion (SVN). However this is not encouraged unless you are
a developer. To check out the SVN source code, you should use the following 
command:

  svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/mmc/trunk mmc

then type the password as "anonymous_user". This will allow you to 
anonymously check out the entire source code tree.

To compile the software, you need to install GNU gcc compiler toolchain
on your system. For Debian/Ubuntu based GNU/Linux systems, you can type

  sudo apt-get install gcc

and for Fedora/Redhat based GNU/Linux systems, you can type

  su -c 'yum install gcc'
 
To compile the binary with multi-threaded computing via OpenMP, 
your gcc version should be at least 4.0. To compile the binary 
supporting SSE4 instructions, gcc version should be at least 
4.3.4. For windows users, you should install MinGW with a later 
version of gcc [3]. You should also install LibGW32C library [4] 
and copy the missing header files from GnuWin32\include\glibc
to MinGW\include when you compile the code (these files typically include
ieee754.h, features.h, endian.h, bits/, gnu/, sys/cdefs.h, sys/ioctl.h 
and sys/ttydefaults.h). For Mac OS X users, you need to install the
mp-gcc4.x series from MacPorts and use the instructions below to compile
the MMC source code.

To compile the program, you should first navigate into the mmc/src folder,
and type

  make release

this will compile a single-threaded, optimized binary under mmc/src/bin
folder. Other "make" options include

  make omp      # this compiles a multi-threaded binary using OpenMP
  make prof     # this makes a binary to produce profiling info for gprof
  make sse      # this uses SSE4 for all vector operations (dot, cross), implies omp
  make ssemath  # this uses SSE4 for both vector operations and math functions
  make          # this produces an non-optimized binary with debugging symbols

If you append "-f makefile_sfmt" at the end of any of the above 
make commands, you will get an executable named "mmc_sfmt", which uses a 
fast MT19937 random-number-generator (RNG) instead of the default GLIBC 
48bit RNG. If your CPU supports SSE4, the fastest binary can be obtained
by running the following command:

  make ssemath -f makefile_sfmt

You should be able to compile the code with an Intel C++ compiler,
an AMD C compiler or LLVM compiler without any difficulty. To use other
compilers, you simply append "CC=compiler_exe" to the above make 
commands. If you see any error messages, please google and fix 
your compiler settings or install the missing libraries.

A special note for Mac OS users: you need to install mp-gcc4{4,5,6}
from MacPorts in order to compile MMC. The default gcc (4.2) installed
by Xcode 3.x does not support thread-local storage. Once downloaded
and installed MacPorts from www.macports.org, you can install gcc by

  sudo port install mp-gcc44

Then add /opt/local/bin to your $PATH variable. A example compilation 
command for MMC looks like

  make ssemath -f makefile_sfmt CC=gcc-mp-4.4

After compilation, you may add the path to the "mmc" binary (typically,
mmc/src/bin) to your search path. To do so, you should modify your 
$PATH environment variable. Detailed instructions can be found at [5].

-------------------------------------------------------------------------------

III. Running Simulations


3.1 Preparation

Before you create/run your own MMC simulations, we suggest you
first understanding all the examples under the mmc/example 
directory, checking out the formats of the input files and the 
scripts for pre- and post-processing.

Because MMC uses FE meshes in the simulation, you should create
a mesh for your problem domain before launching any simulation.
This can be done fairly straightforwardly using a Matlab/Octave 
mesh generator, iso2mesh [1], developed by the MMC author. In 
the mmc/matlab folder, we also provide additional functions to 
generate regular grid-shaped tetrahedral meshes.

It is required to use the "savemmcmesh" function under the 
mmc/matlab folder to save the mesh output from iso2mesh, because 
it performs additional tests to ensure the consistency of element 
orientations. If you choose not to use savemmcmesh, you 
MUST call the "meshreorient" function in iso2mesh to test 
the "elem" array and make sure all elements are oriented in the 
same direction. Otherwise, MMC will give incorrect results.


3.2 Command line options

The full command line options of MMC include the following:
<pre>
usage: mmc <param1> <param2> ...
where possible parameters include (the first item in [] is the default value)
 -i 	       (--interactive) interactive mode
 -s sessionid  (--session)     a string used to tag all output file names
 -f config     (--input)       read config from a file
 -n [0.|float] (--photon)      total photon number, max allowed value is 2^32-1
 -b [0|1]      (--reflect)     1 do reflection at int&ext boundaries, 0 no ref.
 -e [0.|float] (--minenergy)   minimum energy level to trigger Russian roulette
 -U [1|0]      (--normalize)   1 to normalize the fluence to unitary,0 save raw
 -d [0|1]      (--savedet)     1 to save photon info at detectors,0 not to save
 -m [0|1]      (--momentum)    1 to save photon momentum transfer,0 not to save
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save
 -C [1|0]      (--basisorder)  1 piece-wise-linear basis for fluence,0 constant
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh
 -O [X|XFE]    (--outputtype)  X - output flux, F - fluence, E - energy deposit
 -u [1.|float] (--unitinmm)    define the length unit in mm for the mesh
 -h            (--help)        print this message
 -l            (--log)         print messages to a log file instead
 -E [0|int]    (--seed)        set random-number-generator seed
 -M [H|PHBS]   (--method)      choose ray-tracing algorithm (only use 1 letter)
                               P - Plucker-coordinate ray-tracing algorithm
			       H - Havel's SSE4 ray-tracing algorithm
			       B - partial Badouel's method (used by TIM-OS)
			       S - branch-less Badouel's method with SSE
 -D [0|int]    (--debug)       print debug information (you can use an integer
  or                           or a string by combining the following flags)
 -D [''|MCBWDIOXATRPE]         1 M  photon movement info
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
      combine multiple items by using a string, or add selected numbers together
example:
       mmc -n 1000000 -f input.inp -s test -b 0 -D TP
</pre>


3.3 Input files

The simplest example can be found under the "example/onecube" 
folder. Please run "createmesh.m" first from Matlab/Octave to 
create all the mesh files, which include

  elem_onecube.dat    -- tetrahedral element file
  facenb_onecube.dat  -- element neighbors of each face
  node_onecube.dat    -- node coordinates
  prop_onecube.dat    -- optical properties of each element type
  velem_onecube.dat   -- volume of each element

The input file of the example is named "onecube.inp", where we
specify most of the simulation parameters. The input file follows
a similar format as in MCX, which looks like the following

 100                  # total photon number (can be overwriten by -n)
 17182818             # RNG seed, negative to regenerate
 2. 8. 0.             # source position (mm)
 0. 0. 1.             # initial incident vector
 0.e+00 5.e-09 5e-10  # time-gates(s): start, end, step
 onecube              # mesh id: name stub to all mesh files
 3                    # index of element (starting from 1) which encloses the source
 3       1.0          # detector number and radius (mm)
 2.0     6.0    0.0   # detector 1 position (mm)
 2.0     4.0    0.0   # ...
 2.0     2.0    0.0

The mesh files are linked through the "mesh id" (a name stub) with a 
format of {node|elem|facenb|velem}_meshid.dat. All mesh files must 
exist for an MMC simulation. If the index to the tetrahedron that 
encloses the source is not known, please use the "tsearchn" 
function in matlab/octave to find out and supply it in the 7th line
in the input file. Examples are provided in mmc/examples/meshtest/createmesh.m.

To run a simulation, you should execute the "run_test.sh" bash
script in this folder. If you want to run mmc directly from the 
command line, you can do so by typing

 ../../src/bin/mmc -n 20 -f onecube.inp -s onecube 

where -n specifies the total photon number to be simulated,
-f specifies the input file, and -s gives the output file name.
To see all the supported options, run "mmc" without any parameters.

The above command only simulates 20 photons and will complete
instantly. An output file "onecube.dat" will be saved to record the
normalized (unitary) flux at each node. If one specifies
multiple time-windows from the input file, the output will 
contain multiple blocks with each block corresponding to the
time-domain solution at all nodes computed for each time window.

More sophisticated examples can be found under the
example/validation and example/meshtest folders, where you
can find "createmesh" scripts and post-processing script to make
plots from the simulation results.


3.4 JSON-formatted input files

Starting from version 0.9, MMC accepts a JSON-formatted input file in
addition to the conventional tMCimg-like input format. JSON 
(JavaScript Object Notation) is a portable, human-readable and 
"fat-free" text format to represent complex and hierarchical data.
Using the JSON format makes a input file self-explanatory, extensible
and easy-to-interface with other applications (like MATLAB).

A sample JSON input file can be found under the examples/onecube
folder. The same file, onecube.json, is also shown below:

 {
    "Mesh": {
	"MeshID": "onecube",
	"InitElem": 3
    },
    "Session": {
	"Photons":  100,
	"Seed":     17182818,
	"ID":       "onecube"
    },
    "Forward": {
	"T0": 0.0e+00,
	"T1": 5.0e-09,
	"Dt": 5.0e-10
    },
    "Optode": {
	"Source": {
	    "Pos": [2.0, 8.0, 0.0],
	    "Dir": [0.0, 0.0, 1.0]
	},
	"Detector": [
	    {
		"Pos": [2.0, 6.0, 0.0],
		"R": 1.0
	    },
            {
                "Pos": [2.0, 4.0, 0.0],
                "R": 1.0
            },
            {
                "Pos": [2.0, 2.0, 0.0],
                "R": 1.0
            }
	]
    }
 }

A JSON input file requires 4 root objects, namely "Mesh", "Session", "Forward" 
and "Optode". Each object is a data structure providing information
as indicated by its name. Each object can contain various sub-fields. 
The orders of the fields in the same level are flexible. For each field, 
you can always find the equivalent fields in the *.inp input files. 
For example, The "MeshID" field under the "Mesh" object 
is the same as Line#6 in onecube.inp; the "InitElem" under "Mesh" is
the same as Line#7; the "Forward.T0" is the same as the first number 
in Line#5, etc.

An MMC JSON input file must be a valid JSON text file. You can validate
your input file by running a JSON validator, for example http://jsonlint.com/
You should always use "..." to quote a "name" and separate parallel
items by ",".

MMC accepts an alternative form of JSON input, but using it is not 
recommended. In the alternative format, you can use 
 "rootobj_name.field_name": value 
to represent any parameter directly in the root level. For example

 {
    "Mesh.MeshID": "onecube",
    "Session.ID": "onecube",
    ...
 }

You can even mix the alternative format with the standard format. 
If any input parameter has values in both formats in a single input 
file, the standard-formatted value has higher priority.

To invoke the JSON-formatted input file in your simulations, you 
can use the "-f" command line option with MMC, just like using an 
.inp file. For example:

  ../../src/bin/mmc -n 20 -f onecube.json -s onecubejson -D M

The input file must have a ".json" suffix in order for MMC to 
recognize. If the input information is set in both command line,
and input file, the command line value has higher priority
(this is the same for .inp input files). For example, when 
using "-n 20", the value set in "Session"/"Photons" is overwritten 
to 20; when using "-s onecubejson", the "Session"/"ID" value is modified.
If your JSON input file is invalid, MMC will quit and point out
where it expects you to double check.

-------------------------------------------------------------------------------

IV. Plotting the Results

As described above, MMC produces a flux/fluence output file as
"session-id".dat. By default, this file contains the normalized,
i.e. under unitary source, flux at each node of the mesh. The detailed
interpretation of the output data can be found in [6]. If multiple 
time-windows are defined, the output file will contain
multiple blocks of data, with each block being the flux distribution
at each node at the center point of each time-window. The total
number of blocks equals to the total time-gate number.

To read the mesh files (tetrahedral elements and nodes) into matlab, 
one can use readmmcnode and readmmcelem function under the mmc/matlab
directory. Plotting non-structural meshes in matlab is possible with
interpolation functions such as griddata3. However, it is very
time-consuming for large meshes. In iso2mesh, a fast mesh slicing
& plotting function, qmeshcut, is very efficient in making 3D
plots of mesh or cross-sections. More details can be found at 
this webpage [7], or "help qmeshcut" in matlab. Another useful
function is plotmesh in iso2mesh toolbox. It has very flexible
syntax to allow users to plot surfaces, volumetric meshes and
cross-section plots. One can use something like

  plotmesh(node,elem,'x<30 & y>30');

to plot a sliced mesh.

Please edit or browse the *.m files under all example subfolder
to find more options to make plot from MMC output.

When users specify "-d 1" to record partial path lengths for all
detected photons, an output file named "sessionid".mch will be 
saved under the same folder. This file can be loaded into 
Matlab/Octave using the "loadmch.m" script under the mmc/matlab
folder. The output of loadmch script has the following columns:

  detector-id, scattering-events, partial-length_1, partial-length_2, ...., additional data ...

The simulation settings will be returned by a structure. Using the
information from the mch file will allow you to re-scale the detector
readings without rerunning the simulation (for absorption changes only).

-------------------------------------------------------------------------------

V. Known issues and TODOs

* MMC only supports linear tetrahedral elements at this point. Quadratic \
 elements will be added later
* Currently, this code only supports element-based optical properties; \
 nodal-based optical properties (for continuously varying media) will be \
 added in a future release

-------------------------------------------------------------------------------
VI.   Getting Involved

MMC is an open-source software. It is released under the terms of GNU 
General Public License version 3 (GPLv3). That means not only everyone 
can download and use MMC for any purposes, but also you can modify the 
code and share the improved software with others (as long as the derived 
work is also licensed under the GPLv3 license).

If you already made a change to the source code to fix a bug you encountered 
in your research, we are appreciated if you can share your changes (as 
"git diff" outputs) with the developers. We will patch the code as soon 
as we fully test the changes (we will acknowledge your contribution in 
the MMC documentation). If you want to become a developer, please send 
an email to Qianqian and we will review your request. Once permitted, 
you will have developer access to the source code repository.

In you are a user, please use our mmc-users mailing list to post 
questions or share experience regarding MMC. The mailing lists can be
found from this link:

 http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Portal

-------------------------------------------------------------------------------
VII.  Reference

[1] http://iso2mesh.sf.net  -- an image-based surface/volumetric mesh generator
[2] http://mcx.sf.net       -- Monte Carlo eXtreme: a GPU-accelerated MC code
[3] http://sourceforge.net/projects/mingw/files/GCC%20Version%204/
[4] http://developer.apple.com/mac/library/releasenotes/DeveloperTools/RN-llvm-gcc/index.html
[5] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/AddPath
[6] http://mcx.sf.net/cgi-bin/index.cgi?MMC/Doc/FAQ#How_do_I_interpret_MMC_s_output_data
[7] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?fun/qmeshcut
