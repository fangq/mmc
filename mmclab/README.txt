= MMCLAB: MMC for MATLAB and GNU Octave =

*Author: Qianqian Fang <q.fang at neu.edu>
*License: GNU General Public License version 3 (GPLv3)
*Version: this package is part of Mesh-based Monte Carlo (MMC) 3.0, v2026
*URL: https://mcx.space/wiki/?Learn#mmc

<toc>


== # Introduction ==

MMC is a mesh-based Monte Carlo (MC) photon simulation software. It can 
utilize a tetrahedral mesh to model a complex anatomical structure, 
thus, has been shown to be more accurate and computationally efficient 
than the conventional MC codes.

MMCLAB is the native MEX version of MMC for MATLAB and GNU Octave. By
converting the input and output files into convenient in-memory variables, 
MMCLAB is very intuitive to use and straightforward to be integrated with 
mesh generation and post-simulation analyses.

Because MMCLAB contains the same computational codes for multi-threading
photon simulation as in a MMC binary, running MMCLAB inside MATLAB is expected
to give similar speed as running a standalone MMC binary. If your CPU
supports, running mmclab with the 'sse' option can be 25% faster than
the standard mode.


== # Installation ==

Installation of MMCLAB is straightforward. You first download the
MMCLAB package and unzip it to a folder; then you add the folder path 
into MATLAB's search path list. This can be done with the "addpath" 
command in a working session; if you want to add this path permanently,  
use the "pathtool" command, or edit your startup.m (~/.octaverc for 
Octave).

After installation, please type "help mmclab" in MATLAB/Octave to print the
help information.

The Octave mex files released for Windows and macOS require Octave 10 or
newer; the Linux Octave mex file works with any recent Octave version. The
released Linux mex files are built against glibc 2.17 and run on most Linux
distributions released since 2014. MMCLAB can also be installed with the
all-in-one MCX Studio package. For Python users, the same functionality is
provided by pmmc ("pip install pmmc"), whose pmmc.mmclab() function accepts
the same cfg structure (as a dict).


== # How to use MMCLAB in MATLAB/Octave ==

To learn the basic usage of MMCLAB, you can type

  help mmclab

and enter in MATLAB/Octave to see the help information regarding how to use this 
function. The help information is listed below. You can find the input/output 
formats and examples. The input cfg structure has very similar field names as
the verbose command line options in MMC.

<pre>
  #############################################################################%
                      Mesh-based Monte Carlo (MMC) - OpenCL                   %
           Copyright (c) 2010-2026 Qianqian Fang <q.fang at neu.edu>          %
               https://mcx.space/#mmc  &  https://neurojson.io                %
                                                                              %
  Computational Optics & Translational Imaging (COTI) Lab  [https://fanglab.org]%
    Department of Bioengineering, Northeastern University, Boston, MA, USA    %
  #############################################################################%
     The MCX Project is funded by the NIH/NIGMS under grant R01-GM114365      %
  #############################################################################%
   Open-source codes and reusable scientific data are essential for research, %
  MCX proudly developed human-readable JSON-based data formats for easy reuse.%
                                                                              %
  Please visit our free scientific data sharing portal at https://neurojson.io%
  and consider sharing your public datasets in standardized JSON/JData format %
  #############################################################################%
  $Rev::10da53$v2026$Date::2024-03-14 16:10:25 -04$by $Author::Qianqian Fang%
  #############################################################################%
 
  Format:
     [fluence,detphoton,ncfg,seeds,traj]=mmclab(cfg);
           or
     fluence=mmclab(cfg);
     newcfg=mmclab(cfg,'prep');
     [fluence,detphoton,ncfg,seeds,traj]=mmclab(cfg, options);
 
  Input:
     cfg: a struct, or struct array. Each element in cfg defines
          a set of parameters for a simulation.
 
     option: (optional), options is a string, specifying additional options
          option='preview': this plots the domain configuration using mcxpreview(cfg)
          option='opencl':  force using OpenCL (set cfg.gpuid=1 if not set)
                            instead of SSE on CPUs/GPUs that support OpenCL
 
     if one defines USE_MCXCL in the MATLAB/Octave base workspace, mmclab
     inspects it before every call:
         USE_MCXCL = 1   (or unset)   default: run on the first GPU/device
         USE_MCXCL = 0                identical to the default for mmclab
                                       (mmc.mex is unified - no CUDA/OpenCL
                                       switch). Provided for mcxlabcl parity.
         USE_MCXCL = N  with |N| > 1  override cfg.gpuid with abs(N) so that
                                       every cfg in the input runs on device
                                       N (1-based, matches gpuinfo numbering).
                                       Useful for selecting one specific GPU
                                       on multi-GPU systems without editing
                                       every script's cfg.gpuid.
 
 
     cfg may contain the following fields:
 
  == Required ==
      *cfg.nphoton:     the total number of photons to be simulated (integer)
      *cfg.prop:        an N by 4 array, each row specifies [mua, mus, g, n] in order.
                        the first row corresponds to medium type 0 which is
                        typically [0 0 1 1]. The second row is type 1, and so on.
      *cfg.node:        node array for the input tetrahedral mesh, 3 columns: (x,y,z)
      *cfg.elem:        element array for the input tetrahedral mesh, 4 columns
      *cfg.elemprop:    element property index for input tetrahedral mesh
      cfg.face:        (OptiX only) a surface-only mesh used instead of cfg.elem, 5 columns:
                       3 node indices (1-based), the inner and the outer medium labels;
                       the triangle normal (right-hand rule) points from inner to outer
      *cfg.tstart:      starting time of the simulation (in seconds)
      *cfg.tstep:       time-gate width of the simulation (in seconds)
      *cfg.tend:        ending time of the simulation (in second)
      *cfg.srcpos:      a 1 by 3 vector, the position of the source in mesh node length unit
      *cfg.srcdir:      if defined as [vx, vy, vy], it specifies the incident vector
                        if defined as [vx, vy, vy, focus], the first 3 elements define
                        the incident vector; focus controls the convergence or
                        divergence of the beam:
                        focus=0: collimated beam
                        focus<0: diverging beam from an imaginary src at c0-|focus|*[vx vy vz]
                        focus>0: converging beam, focusing to a point at c0+|focus|*[vx vy vz]
                        where c0 is the centroid of the source domain. Setting focus does
                        not impact pencil/isotropic/cone sources.
 
  == MC simulation settings ==
       cfg.seed:        seed for the random number generator (integer)
                        if set to a uint8 array, the binary data in each column is used
                        to seed a photon (i.e. the "replay" mode), default value: 1648335518
       cfg.isreflect:   [1]-consider refractive index mismatch, 0-matched index
                        2 - total absorption on exterior surface
                        3 - prefect reflection (mirror) on exterior surface
       cfg.isnormalized:[1]-normalize the output fluence to unitary source, 0-no reflection
       cfg.isspecular:  [1]-calculate specular reflection if source is outside
       cfg.ismomentum:  [0]-save momentum transfer for each detected photon
       cfg.method:      ray-tracing method, ["plucker"]:Plucker, "havel": Havel (SSE4),
                        "badouel": partial Badouel, "elem": branchless Badouel (SSE),
                        "grid": dual-grid MMC
       cfg.mcmethod:    0 use MCX-styled MC method, 1 use MCML style MC
       cfg.nout:        [1.0] refractive index for medium type 0 (background)
       cfg.minenergy:   terminate photon when weight less than this level (float) [0.0]
       cfg.roulettesize:[10] size of Russian roulette
       cfg.steps:       [dx, dy, dz], defines the DMMC grid voxel size,
                        must be isotropic, i.e. dx=dy=dz, only used when
                        cfg.method = 'grid', by default dx=dy=dz=1
       cfg.unitinmm:    defines the default length unit (to interpret mesh nodes, src/det positions
                        the default value is 1.0 (mm). For example, if the mesh node length unit is
                        in cm, one should set unitinmm to 10.
       cfg.basisorder:  [1]-linear basis, 0-piece-wise constant basis
 
  == Source-detector parameters ==
       cfg.detpos:      an N by 4 array, each row specifying a detector: [x,y,z,radius]
       cfg.detparam1:   a 4-element vector: [u1, u2, u3, Nx], define the x-axis of a widefield
                        detector, [u1,u2,u3] is the x-axis vector; Nx is the pixel count in the x-axis
       cfg.detparam2:   a 4-element vector: [v1, v2, v3, Ny], define the y-axis of a widefield
                        detector, [v1,v2,v3] is the y-axis vector; Ny is the pixel count in the y-axis
       cfg.maxdetphoton:   maximum number of photons saved by the detectors [1000000]
       cfg.srctype:     source type, the parameters of the src are specified by cfg.srcparam{1,2}
                       'pencil' - default, pencil beam, no param needed
                       'isotropic' - isotropic source, no param needed
                       'cone' - uniform cone beam, srcparam1(1) is the half-angle in radian
                       'gaussian' - a gaussian beam, srcparam1(1) specifies the waist radius
                                 (in default length unit); if one specifies a non-zero focal length
                                 using cfg.srcdir, the gaussian beam can be converging to or
                                 diverging from the waist center, which is located at srcpos+focus*srcdir;
                                 optionally, one can specify the wavelength lambda (in cfg.unitinmm mm),
                                 using srcparam1(2). This will rescale the Gaussian profile according
                                 to w(z)=w0*sqrt(1-(z/z0)^2), where w0 is the waist radius, z is the
                                 distance (in mm) to the waist center (focus), and z0 is the Rayleigh
                                 range (in mm), and z0 is related to w0 by z0=w0^2*pi/lambda
                       'planar' - a 3D quadrilateral uniform planar source, with three corners specified
                                 by srcpos, srcpos+srcparam1(1:3) and srcpos+srcparam2(1:3)
                       'pattern' - a 3D quadrilateral pattern illumination, same as above, except
                                 srcparam1(4) and srcparam2(4) specify the pattern array x/y dimensions,
                                 and srcpattern is a double-precision pattern array, with values between [0-1].
                                 if cfg.srcnum>1, srcpattern must be a 3-D double-precision array with
                                 a dimension of [srcnum srcparam1(4) srcparam2(4)]
                                 Example: <demo_photon_sharing.m>
                       'fourier' - spatial frequency domain source, similar to 'planar', except
                                 the integer parts of srcparam1(4) and srcparam2(4) represent
                                 the x/y frequencies; the fraction part of srcparam1(4) multiplies
                                 2*pi represents the phase shift (phi0); 1.0 minus the fraction part of
                                 srcparam2(4) is the modulation depth (M). Put in equations:
                                     S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
                       'arcsine' - similar to isotropic, except the zenith angle is uniform
                                 distribution, rather than a sine distribution.
                       'disk' - a uniform disk source pointing along srcdir; the radius is
                                set by srcparam1(1) (in default length unit)
                       'fourierx' - a general Fourier source, the parameters are
                                srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phi0,M]
                                normalized vectors satisfy: srcdir cross v1=v2
                                the phase shift is phi0*2*pi
                       'fourierx2d' - a general 2D Fourier basis, parameters
                                srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phix,phiy]
                                the phase shift is phi{x,y}*2*pi
                       'zgaussian' - an angular gaussian beam, srcparam1(1) specifies the variance in
                                the zenith angle
       cfg.{srcparam1,srcparam2}: 1x4 vectors, see cfg.srctype for details
       cfg.srcpattern: see cfg.srctype for details
       cfg.srcnum:     the number of source patterns that are
                       simultaneously simulated; only works for 'pattern'
                       source, see cfg.srctype='pattern' for details
                       Example <demo_photon_sharing.m>
       cfg.replaydet:  only works when cfg.outputtype is 'jacobian', 'wl', 'nscat', or 'wp' and cfg.seed is an array
                        0 replay all detectors and sum all Jacobians into one volume
                        a positive number: the index of the detector to replay and obtain Jacobians
       cfg.voidtime:   for wide-field sources, [1]-start timer at launch, 0-when entering
                       the first non-zero voxel
 
       by default, mmc assumes the mesh and source position settings are all in mm unit.
       if the mesh coordinates/source positions are not in mm unit, one needs to define
       cfg.unitinmm  (in mm) to specify the actual length unit.
 
  == Optional mesh data ==
      -cfg.facenb:      element face neighborhood list (calculated by faceneighbors())
      -cfg.evol:        element volume (calculated by elemvolume() with iso2mesh)
      -cfg.e0:          the element ID enclosing the source, if not defined,
                        it will be calculated by tsearchn(node,elem,srcpos);
                        if cfg.e0 is set as one of the following characters,
                        mmclab will do an initial ray-tracing and move
                        srcpos to the first intersection to the surface:
                        '>': search along the forward (srcdir) direction
                        '<': search along the backward direction
                        '-': search both directions
 
  == Output control ==
       cfg.issaveexit: [0]-whether to save the position (x,y,z) and (vx,vy,vz) for a detected photon
                       When issaveexit is set to 2, this signals mmclab to output a 2D array
                       a widefield detector defined by cfg.detparam1 and cfg.detparam2.
                       Specifically, instead of output individual detected photon data,
                       the exiting weight of all detected photons are binned and summed
                       into a detparam1(4)-by-detparam2(4) pixelated output, covering
                       the quadrilateral area with 4 corners defined by (only 1:3 elem are used)
                       [detpos, detpos+detparam1, detpos+detparam1+detparam2, detpos+detparam2].
                       When a widefield detector is used, detector radii are ignored.
       cfg.issaveref:  [0]-save diffuse reflectance/transmittance on the exterior surfaces.
                       The output is stored as flux.dref in a 2D array of size [#Nf,  #time_gate]
                       where #Nf is the number of triangles on the surface; #time_gate is the
                       number of total time gates. To plot the surface diffuse reflectance, the output
                       triangle surface mesh can be extracted by faces=faceneighbors(cfg.elem,'rowmajor');
                       where 'faceneighbors' can be found in the iso2mesh toolbox.
                       Example: see <demo_mmclab_basic.m>
       cfg.issaveseed:  [0]-save the RNG seed for a detected photon so one can replay
       cfg.isatomic:    deprecated and ignored: fluence is always saved with atomic operations
       cfg.outputtype:  'flux' - output fluence-rate
                        'fluence' - fluence,
                        'energy' - energy deposit,
                        'jacobian' - mua Jacobian (replay mode)
                        'wl'- weighted path lengths to build mua Jacobian (replay mode)
                        'wp'- weighted scattering counts to build mus Jacobian (replay mode)
       cfg.debuglevel:  debug flag string, a subset of [SCBWDIOXATRPEM], no space
       cfg.debugphoton: print the photon movement debug info only for a specified photon ID
       cfg.maxjumpdebug: [10000000|int] when trajectory is requested in the output,
                      use this parameter to set the maximum position stored. By default,
                      only the first 1e6 positions are stored.
 
       fields marked with * are required; options in [] are the default values
       fields marked with - are calculated if not given (can be faster if precomputed)
 
 
     type: omit or 'omp' for multi-threading version; 'sse' for the SSE4 MMC,
           the SSE4 version is about 25% faster, but requires newer CPUs;
           if type='prep' with a single output, mmclab returns ncfg only.
 
  Output:
       fluence: a struct array, with a length equals to that of cfg.
             For each element of fluence, fluence(i).data is a 2D array with
             dimensions [size(cfg.node,1), total-time-gates] if cfg.basisorder=1,
             or [size(cfg.elem,1), total-time-gates] if cfg.basisorder=0.
             The content of the array is the normalized fluence-rate (or others
             depending on cfg.outputtype) at each mesh node and time-gate.
 
             If cfg.issaveref is set to 1, fluence(i).dref is not empty, and stores
             the surface diffuse reflectance (normalized by default). The surface mesh
             that the dref output is attached can be obtained by faces=faceneighbors(cfg.elem,'rowmajor');
       detphoton: (optional) a struct array, with a length equals to that of cfg.
             Starting from v2016.5, the detphoton contains the below subfields:
               detphoton.detid: the ID(>0) of the detector that captures the photon
               detphoton.nscat: cumulative scattering event counts in each medium
               detphoton.ppath: cumulative path lengths in each medium (partial pathlength)
                    one need to multiply cfg.unitinmm with ppath to convert it to mm.
               detphoton.mom: cumulative cos_theta for momentum transfer in each medium
               detphoton.p or .v: exit position and direction, when cfg.issaveexit=1
               detphoton.w0: photon initial weight at launch time
               detphoton.prop: optical properties, a copy of cfg.prop
               detphoton.data: a concatenated and transposed array in the order of
                     [detid nscat ppath mom p v w0]'
               "data" is the is the only subfield in all MMCLAB before 2016.5
       ncfg: (optional), if given, mmclab returns the preprocessed cfg structure,
             including the calculated subfields (marked by "-"). This can be
             used in the subsequent simulations to avoid repetitive preprocessing.
       seeds: (optional), if give, mmclab returns the seeds, in the form of
             a byte array (uint8) for each detected photon. The column number
             of seed equals that of detphoton.
       trajectory: (optional), if given, mmclab returns the trajectory data for
             each simulated photon. The output has 6 rows, the meanings are
                id:  1:    index of the photon packet
                pos: 2-4:  x/y/z/ of each trajectory position
                     5:    current photon packet weight
                     6:    enclosing element's ID
             By default, mcxlab only records the first 1e7 positions along all
             simulated photons; change cfg.maxjumpdebug to define a different limit.
 
  Example:
       cfg.nphoton=1e5;
       [cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
       cfg.elemprop=ones(size(cfg.elem,1),1);
       cfg.srcpos=[30 30 0];
       cfg.srcdir=[0 0 1];
       cfg.prop=[0 0 1 1;0.005 1 0 1.37];
       cfg.tstart=0;
       cfg.tend=5e-9;
       cfg.tstep=5e-10;
       cfg.debuglevel='TP';
       % preprocessing to populate the missing fields to save computation
       ncfg=mmclab(cfg,'prep');
 
       cfgs(1)=ncfg;   % when using struct array input, all fields must be defined
       cfgs(2)=ncfg;
       cfgs(1).isreflect=0;
       cfgs(2).isreflect=1;
       cfgs(2).detpos=[30 20 0 1;30 40 0 1;20 30 1 1;40 30 0 1];
       % calculate the fluence and partial path lengths for the two configurations
       [fluxs,detps]=mmclab(cfgs);
 
 
  This function is part of Mesh-based Monte Carlo (MMC) URL: https://mcx.space/#mmc
 
  License: GNU General Public License version 3, please read LICENSE.txt for details
</pre>

== # Examples ==

We provided several examples to demonstrate the basic usage of MMCLAB;
some of the examples also perform validations of MMC algorithm in both
homogeneous and heterogeneous domains. These examples are explained below:

==== demo_mmclab_basic.m ====

In this example, we show the most basic usage of MMCLAB. This includes
how to define the input configuration structure, launch MMC simulations
and plot the output data.

==== demo_example_validation.m ====

In this example, we validate MMCLAB with a homogeneous medium in a 
cubic domain. This is the same as in mmc/examples/validation;
the details of the simulation is described in Fig.2 of [Fang2010].

==== demo_example_meshtest.m ====

In this example, we validate the MMCLAB solver with a heterogeneous
domain and the analytical solution from the diffusion model. The 
domain is consisted of a 6x6x6 cm box with a 2cm diameter sphere 
embedded at the center. 

This test is identical to the simulations under mmc/examples/meshtest, 
which are described by Fig. 3 in [Fang2010].

==== demo_example_onecube.m ====

In this example, we demonstrate how to use the debug flags
to print photon trajectories. The domain is a simple 10x10x10mm
cube.

This test is similar to the simulations under mmc/examples/onecube. 

==== demo_example_replay.m ====

In this example, we run MMC using a homogeneous cubic domain,
same as in demo_example_validation.m. We save the seeds of the
detected photons, and then rerun these photons again (i.e.
the "replay" part). This feature allows one to identify features
that are only related to a source/detector pair.

==== demo_compare_mmc_mcx.m ====

In this example, we compare MMC and MCX using both simple and
widefield sources.

==== demo_sfdi_2layer.m ====

In this example, we simulate an spatial-frequency domain imaging source
using a 2-layer brain model.

==== demo_mmclab_slit.m ====

In this example, we show how to use a slit (line light source) in
an MMC simulation.


== # How to compile MMCLAB ==

To compile MMCLAB from source code, you need to make sure your 
computer have the following requirements:

* your computer should have the MATLAB C compiler (mex) and/or the Octave development package (liboctave-dev on Ubuntu) installed
* to compile the SSE versions, your computer should support SSE4 instructions
* the commands "mex" and "mkoctfile" should be in the search path ($PATH)
* you should have both gcc and g++ with version 4.4 or newer
* for windows, you should install Cygwin64 with mingw64-i686-gcc and make packages
* for windows, you need to create a mexopts.bat file as \
C:\Users\<username>\AppData\Roaming\MathWorks\MATLAB\RXXXX\mexopts.bat \
with content modified from mmc/src/mexopts_cygwin64_gcc.bat

To compile MMCLAB for MATLAB, you need to cd mmc/src directory, and type 

 make mex
or
 make mexsse 

from a shell window. You need to make sure your MATLAB is installed and 
the command <tt>mex</tt> is included in your PATH environment variable. Similarly, 
to compile MMCLAB for Octave, you type

 make oct
or
 make octsse

The command <tt>mkoctfile</tt> must be accessible from your command line;
on Ubuntu/Debian, it is provided by the package <tt>liboctave-dev</tt>.

Compiling MMCLAB on Windows and Mac OS, you should use the following
command:

 make ... -f makefile_sfmt ...

If you have multiple gcc/g++ installed on your system, and the default 
gcc and g++ is older than 4.4 (use gcc -v to print version), you can 
use the following command to specify the correct gcc command:

 make ... CC=/path/to/gcc/4.4 CXX=/path/to/g++/4.4


== # Screenshots ==

Screenshot for using MMCLAB in MATLAB:
  http://mcx.sourceforge.net/upload/matlab_mmclab.png

Screenshot for using MMCLAB in GNU Octave:
  http://mcx.sourceforge.net/upload/octave_mmclab.png


== # Reference ==

[Fang2019] Qianqian Fang* and Shijie Yan, "GPU-accelerated mesh-based \
 Monte Carlo photon transport simulations," J. of Biomedical Optics, in press, 2019. \
 Preprint URL: https://www.biorxiv.org/content/10.1101/815977v1


[Brain2Mesh2020] Anh Phong Tran† , Shijie Yan† , Qianqian Fang*, (2020) "Improving \
 model-based fNIRS analysis using mesh-based anatomical and light-transport models," \ 
 Neurophotonics, 7(1), 015008, URL: https://doi.org/10.1117/1.NPh.7.1.015008

[Fang2010] Fang Q, "Mesh-based Monte Carlo method using fast ray-tracing \
   in Plucker coordinates," Biomed. Opt. Express 1(1), 165-175 (2010) 

