===============================================================================
=                       Mesh-based Monte Carlo (MMC)                          =
=                     Multi-threaded Edition with SSE4                        =
===============================================================================

Author:  Qianqian Fang <q.fang at neu.edu>
License: GNU General Public License version 3 (GPL v3), see License.txt
Version: 1.4.8-2 (v2019.4, Pork Rinds - beta, update 2)
URL:     http://mcx.space/mmc

-------------------------------------------------------------------------------

Table of Content:

O.   What's New
I.   Introduction
II.  Downloading and Compiling MMC
III. Running Simulations
IV.  Plotting the Results
V.   Known Issues and TODOs
VI   Getting Involved
VII. Acknowledgement
VIII.Reference


---------------------------------------------------------------------

O.    What's New

In MMC v2019.4 (1.4.8-2), the follow feature was added

* Support -X/--saveref to save diffuse reflectance/transmittance on mesh surface
* Speed up DMMC memory operations

It also fixed the below critical bugs:

* fix #35 - incorect mch file header in photon-sharing implementation
* restore the capability to save mch files without needing --saveexit 1 
* for Win64, use a newer version of libgomp-1.dll to run mmclab without dependency errors


Also, in MMC v2019.3 (1.4.8), we added a list of major new additions, including

* Add 2 built-in complex domain examples - USC_19-5 brain atlas and mcxyz skin-vessel benchmark
* Initial support of "photon sharing" - a fast approach to simultaneouly simulate multiple pattern src/det, as detailed in our Photoncs West 2019 talk by Ruoyang Yao/Shijie Yan [Yao&Yan2019]
* Dual-grid MMC (DMMC) paper published [Yan2019], enabled by "-M G" or cfg.method='grid'
* Add clang compiler support, nightly build compilation script, colored command line output, and more

In addition, we also fixed a number of critical bugs, such as

* fix mmclab gpuinfo output crash using multiple GPUs
* disable linking to Intel OMP library (libiomp5) to avoid MATLAB 2016-2017 crash
* fix a bug for doubling thread number every call to mmc, thanks to Shijie
* fix mmclab crash due to photo sharing update

'''[Yan2019]''' Shijie Yan, Anh Phong Tran, Qianqian Fang*, "A dual-grid mesh-based\
Monte Carlo algorithm for efficient photon transport simulations in complex 3-D media,"\
J. of Biomedical Optics, 24(2), 020503 (2019). URL: https://doi.org/10.1117/1.JBO.24.2.020503

'''[Yao&Yan2019]''' Ruoyang Yao, Shijie Yan, Xavier Intes, Qianqian Fang,  \
"Accelerating Monte Carlo forward model with structured light illumination via 'photon sharing'," \
Photonics West 2019, paper#10874-11, San Francisco, CA, USA. \
[https://www.spiedigitallibrary.org/conference-presentations/10874/108740B/Accelerating-Monte-Carlo-forward-model-with-structured-light-illumination-via/10.1117/12.2510291?SSO=1 Full presentation for our invited talk]

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
  URL: https://www.osapublishing.org/boe/abstract.cfm?uri=boe-1-1-165

While the original MMC paper was based on the Plücker coordinates, a number
of more efficient SIMD-based ray-tracers, namely, Havel SSE4 ray-tracer, 
Badouel SSE ray-tracer and branchless-Badouel SSE ray-tracer (fastest) have 
been added since 2011. These methods can be selected by the -M flag. The 
details of these methods can be found in the below paper

  Qianqian Fang and David R. Kaeli, 
  "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"
  Biomed. Opt. Express 3(12), 3223-3230 (2012)
  URL: https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223
  
and their key differences compared to another mesh-based MC simulator, 
TIM-OS, are discussed in 

  Qianqian Fang, "Comment on 'A study on tetrahedron-based inhomogeneous 
  Monte-Carlo optical simulation'," Biomed. Opt. Express, vol. 2(5) 1258-1264, 2011.
  URL: https://www.osapublishing.org/boe/abstract.cfm?uri=boe-2-5-1258

In addition, the generalized MMC algorithm for wide-field sources and detectors are
described in the following paper, and was made possible with the collaboration
with Ruoyang Yao and Prof. Xavier Intes from RPI

  Yao R, Intes X, Fang Q, "Generalized mesh-based Monte Carlo for
  wide-field illumination and detection via mesh retessellation,"
  Biomed. Optics Express, 7(1), 171-184 (2016)
  URL: https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171

In addition, we have been developing a fast approach to build the
Jacobian matrix for solving inverse problems. The technique is called
"photon replay", and is described in details in the below paper:

  Yao R, Intes X, Fang Q, "A direct approach to compute Jacobians for 
  diffuse optical tomography using perturbation Monte Carlo-based 
  photon 'replay'," Biomed. Optics Express, in press, (2018)

In 2019, we published an improved MMC algorithm, named "dual-grid MMC", 
or DMMC, in the below JBO Letter. This method allows to use separate mesh
for ray-tracing and fluence storage, and can be 2 to 3 fold faster
than the original MMC without loss of accuracy. 

  Shijie Yan, Anh Phong Tran, Qianqian Fang*, "A dual-grid mesh-based 
  Monte Carlo algorithm for efficient photon1transport simulations in 
  complex 3-D media," J. of Biomedical Optics, 24(2), 020503 (2019).

The authors of the papers are greatly appreciated if you can cite 
the above papers as references if you use MMC and related software
in your publication.

-------------------------------------------------------------------------------

II. Download and Compile MMC

The latest release of MMC can be downloaded from the following URL:

  http://mcx.space/#mmc

The development branch (not fully tested) of the code can be accessed 
using Git. However this is not encouraged unless you are
a developer. To check out the Git source code, you should use the following 
command:

  git clone https://github.com/fangq/mmc.git mmc

To compile the software, you need to install GNU gcc compiler toolchain
on your system. For Debian/Ubuntu based GNU/Linux systems, you can type

  sudo apt-get install gcc

and for Fedora/Redhat based GNU/Linux systems, you can type

  su -c 'yum install gcc'
 
To compile the binary with multi-threaded computing via OpenMP, 
your gcc version should be at least 4.0. To compile the binary 
supporting SSE4 instructions, gcc version should be at least 
4.3.4. For windows users, you should install Cygwin64 [3]. During the 
installation, please select mingw64-x86_64-gcc and make packages.
You should also install LibGW32C library [4] and copy the missing 
header files from GnuWin32\include\glibc to C:\cygwin64\usr\include
when you compile the code (these files typically include
ieee754.h, features.h, endian.h, bits/, gnu/, sys/cdefs.h, sys/ioctl.h 
and sys/ttydefaults.h). For Mac OS X users, you need to install the
mp-gcc4.x series from MacPorts or Homebrew and use the instructions 
below to compile the MMC source code.

To compile the program, you should first navigate into the mmc/src folder,
and type

  make

this will create a fully optimized, multi-threaded and SSE4 enabled 
mmc executable, located under the mmc/src/bin/ folder.

Other compilation options include

  make omp      # this compiles a multi-threaded binary using OpenMP
  make release  # create a single-threaded optimized binary
  make prof     # this makes a binary to produce profiling info for gprof
  make sse      # this uses SSE4 for all vector operations (dot, cross), implies omp
  make ssemath  # this uses SSE4 for both vector operations and math functions

if you want to generate a portable binary that does not require external 
library files, you may use (only works for Linux and Windows with gcc)

  make EXTRALIB="-static -lm" # similar to "make", except the binary includes all libraries

if you wish to build the mmc mex file to be used in matlab, you should run

  make mex      # this produces mmc.mex* under mmc/mmclab/ folder

similarly, if you wish to build the mex file for GNU Octave, you should run

  make oct      # this produces mmc.mex* under mmc/mmclab/ folder

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

  make ssemath CC=gcc-mp-4.4

After compilation, you may add the path to the "mmc" binary (typically,
mmc/src/bin) to your search path. To do so, you should modify your 
$PATH environment variable. Detailed instructions can be found at [5].

You can also compile MMC using Intel's C++ compiler - icc. To do this, you run

  make CC=icc

you must enable icc related environment variables by source the compilervars.sh 
file. The speed of icc-generated mmc binary is generally faster than those compiled by 
gcc.

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
###############################################################################
#                         Mesh-based Monte Carlo (MMC)                        #
#          Copyright (c) 2010-2019 Qianqian Fang <q.fang at neu.edu>          #
#                            http://mcx.space/#mmc                            #
#                                                                             #
#Computational Optics & Translational Imaging (COTI) Lab  [http://fanglab.org]#
#            Department of Bioengineering, Northeastern University            #
#                                                                             #
#                Research funded by NIH/NIGMS grant R01-GM114365              #
###############################################################################
$Rev::8270b9$2019.4 $Date::2019-04-24 14:18:58 -04$ by $Author::Qianqian Fang $
###############################################################################

usage: mmc <param1> <param2> ...
where possible parameters include (the first item in [] is the default value)

== Required option ==
 -f config     (--input)       read an input file in .inp or .json format

== MC options ==
 -n [0.|float] (--photon)      total photon number, max allowed value is 2^32-1
 -b [0|1]      (--reflect)     1 do reflection at int&ext boundaries, 0 no ref.
 -U [1|0]      (--normalize)   1 to normalize the fluence to unitary,0 save raw
 -m [0|1]      (--mc)          0 use MCX-styled MC method, 1 use MCML style MC
 -C [1|0]      (--basisorder)  1 piece-wise-linear basis for fluence,0 constant
 -u [1.|float] (--unitinmm)    define the mesh data length unit in mm
 -E [1648335518|int|mch](--seed) set random-number-generator seed;
                               if an mch file is followed, MMC "replays" 
                               the detected photons; the replay mode can be used
                               to calculate the mua/mus Jacobian matrices
 -P [0|int]    (--replaydet)   replay only the detected photons from a given 
                               detector (det ID starts from 1), use with -E 
 -M [H|PHBSG] (--method)      choose ray-tracing algorithm (only use 1 letter)
                               P - Plucker-coordinate ray-tracing algorithm
			       H - Havel's SSE4 ray-tracing algorithm
			       B - partial Badouel's method (used by TIM-OS)
			       S - branch-less Badouel's method with SSE
			       G - dual-grid MMC (DMMC) with voxel data output
 -e [1e-6|float](--minenergy)  minimum energy level to trigger Russian roulette
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh
 -k [1|0]      (--voidtime)    when src is outside, 1 enables timer inside void
 --atomic [1|0]                1 use atomic operations, 0 use non-atomic ones

== Output options ==
 -s sessionid  (--session)     a string used to tag all output file names
 -O [X|XFEJLP] (--outputtype)  X - output flux, F - fluence, E - energy deposit
                               J - Jacobian, L - weighted path length, P -
                               weighted scattering count (J,L,P: replay mode)
 -d [0|1]      (--savedet)     1 to save photon info at detectors,0 not to save
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save
 -x [0|1]      (--saveexit)    1 to save photon exit positions and directions
                               setting -x to 1 also implies setting '-d' to 1
 -X [0|1]      (--saveref)     save diffuse reflectance/transmittance on the 
                               exterior surface. The output is stored in a 
                               file named *_dref.dat, and the 2nd column of 
			       the data is resized to [#Nf, #time_gate] where
			       #Nf is the number of triangles on the surface; 
			       #time_gate is the number of total time gates. 
			       To plot the surface diffuse reflectance, the 
			       output triangle surface mesh can be extracted
			       by faces=faceneighbors(cfg.elem,'rowmajor');
                               where 'faceneighbors' is part of Iso2Mesh.
 -q [0|1]      (--saveseed)    1 save RNG seeds of detected photons for replay
 -F format     (--outputformat)'ascii', 'bin' (in 'double'), 'mc2' (double) 
                               'hdr' (Analyze) or 'nii' (nifti, double)

== User IO options ==
 -h            (--help)        print this message
 -v            (--version)     print MMC version information
 -l            (--log)         print messages to a log file instead
 -i 	       (--interactive) interactive mode

== Debug options ==
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
 --debugphoton [-1|int]        to print the debug info specified by -D only for
                               a single photon, followed by its index (start 0)

== Additional options ==
 --momentum     [0|1]          1 to save photon momentum transfer,0 not to save

== Example ==
       mmc -n 1000000 -f input.json -s test -b 0 -D TP
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
 pencil               # optional: source type
 0 0 0 0              # optional: source parameter set 1
 0 0 0 0              # optional: source parameter set 2

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
    "Domain": {
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
            "Type": "pencil",
	    "Pos": [2.0, 8.0, 0.0],
	    "Dir": [0.0, 0.0, 1.0],
            "Param1": [0.0, 0.0, 0.0, 0.0],
            "Param2": [0.0, 0.0, 0.0, 0.0]
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

A JSON input file requires 4 root objects, namely "Domain", "Session", "Forward" 
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
    "Domain.MeshID": "onecube",
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


3.5 Photon debugging information using -D flag

the output format for -D M (photon moving) is below:

? px py pz eid id scat

? is a single letter representing the state of the current position:
   B a boundary point
   P the photon is passing an interface point
   T the photon terminates at this location due to
      exceeding end of the time window
   M a position other than any of the above

px,py,pz: the current photon position

eid: the index (starting from 1) of the current enclosing element

id: the index of the current photon, from 1 to nphoton

scat: the "normalized" length to read the next scattering site, \
   it is unitless

for -D A (flux accumulation debugging), the output is

A ax ay az ww eid dlen

ax ay az: the location where the accumulation calculation was done \
   (typically, the half-way point of the line segment between the last \
   and current positions)

ww: the photon weight loss for the line segment

dlen=scat/mus of the current element: the distance left to arrive \
   the next scattering site

for -D E

E  px py pz vx vy vz w eid

vx vy vz: the unitary propagation vector when the photon exits
w: the current photon weight

-------------------------------------------------------------------------------

IV. Plotting the Results

As described above, MMC produces a fluence-rate output file as
"session-id".dat. By default, this file contains the normalized,
i.e. under unitary source, fluence at each node of the mesh. The detailed
interpretation of the output data can be found in [6]. If multiple 
time-windows are defined, the output file will contain
multiple blocks of data, with each block being the fluence distribution
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

  plotmesh([node fluence],elem,'x<30 & y>30');

to plot a sliced mesh, or 

  plotmesh([node log10(fluence)],elem,'x=30'); view(3)

to show a cross-sectional plot.

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
the MMC documentation). 

When making edits to the source code with an intent of sharing with the
upstream authors, please set your editor's tab width to 8 so that the 
indentation of the source is correctly displayed. Please keep your patch
as small and local as possible, so that other parts of the code are not
influenced.

To streamline the process process, the best way to contribute your patch
is to click the "fork" button from http://github.com/fangq/mmc, and 
then change the code in your forked repository. Once fully tested and
documented, you can then create a "pull request" so that the upstream author
can review the changes and accept your change.

In you are a user, please use our mmc-users mailing list to post 
questions or share experience regarding MMC. The mailing lists can be
found from this link:

 http://mcx.space/#about

-------------------------------------------------------------------------------

VII.  Acknowledgement

MMC uses the following open-source libraries:

=== SSE Math library by Julien Pommier ===

  Copyright (C) 2007  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)

=== cJSON library by Dave Gamble ===

  Copyright (c) 2009 Dave Gamble

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.

=== SFMT library by Mutsuo Saito, Makoto Matsumoto and Hiroshima University  ===

  Copyright (c) 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
  University. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

      * Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above
	copyright notice, this list of conditions and the following
	disclaimer in the documentation and/or other materials provided
	with the distribution.
      * Neither the name of the Hiroshima University nor the names of
	its contributors may be used to endorse or promote products
	derived from this software without specific prior written
	permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=== drand48_r port for libgw32c by Free Software Foundation ===

   Copyright (C) 1995, 1997, 2001 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, August 1995.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.

=== git-rcs-keywords by Martin Turon (turon) at Github ===

   MMC includes a pair of git filters (.git_filters/rcs-keywords.clean
   and .git_filters/rcs-keywords.smudge) to automatically update SVN
   keywords in mcx_utils.c. The two simple filter scripts were licensed
   under the BSD license according to this link:

   https://github.com/turon/git-rcs-keywords/issues/4

   Both filter files were significantly modified by Qianqian Fang.
 
-------------------------------------------------------------------------------

VIII.  Reference

[1] http://iso2mesh.sf.net  -- an image-based surface/volumetric mesh generator
[2] http://mcx.sf.net       -- Monte Carlo eXtreme: a GPU-accelerated MC code
[3] https://cygwin.com/setup-x86_64.exe
[4] http://developer.apple.com/mac/library/releasenotes/DeveloperTools/RN-llvm-gcc/index.html
[5] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/AddPath
[6] http://mcx.sf.net/cgi-bin/index.cgi?MMC/Doc/FAQ#How_do_I_interpret_MMC_s_output_data
[7] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?fun/qmeshcut
