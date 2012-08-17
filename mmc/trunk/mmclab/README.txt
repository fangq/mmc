= MMCLAB: MMC for MATLAB and GNU Octave =

Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
License: GNU General Public License version 3 (GPLv3)
Version: this package is part of Mesh-based Monte Carlo (MMC) 0.7.9

<toc>


== # Introduction ==

MMC is a mesh-based Monte Carlo (MC) photon simulation software. It uses 
a tetrahedral mesh to model complex anatomical structures, thus, is considered
to be more accurate and computationally efficient than the conventional MC 
codes.

MMCLAB is the native MEX version of MMC for MATLAB and GNU Octave. By
converting the input and output files into convenient in-memory variables, 
MMCLAB is very intuitive to use and straightforward to be integrated with 
mesh-generation and post-simulation analyses.

Because MMCLAB contains the same computational codes for multi-threading
photon simulation as in a MMC binary, running MMCLAB inside MATLAB is expected
to give similar speed as running a standalone MMC binary. 


== # Installation ==

Installation of MMCLAB is very straightforward. You first download and unzip 
the MMCLAB package into a folder; then you add the folder path into MATLAB's
search path list. This can be done with the "addpath" command in a working 
session; if you want to add this path permanently,  use the "pathtool" 
command, or edit your startup.m (~/.octaverc for Octave).

After installation, please type "help mmclab" in MATLAB/Octave to print the
help information.


== # How to use MMCLAB in MATLAB/Octave ==

To learn the basic usage of MMCLAB, you can type

  help mmclab

and enter in MATLAB/Octave to see the help information regarding how to use this 
function. The help information is listed below. You can find the input/output 
formats and examples. The input cfg structure has very similar field names as
the verbose command line options in MMC.

<pre>====================================================================
      MMCLAB - Mesh-based Monte Carlo (MMC) for MATLAB/GNU Octave
--------------------------------------------------------------------
Copyright (c) 2012 Qianqian Fang <fangq at nmr.mgh.harvard.edu>
		   URL: http://mcx.sf.net/mmc/
====================================================================

 Format:
    [flux,detphoton,ncfg]=mmclab(cfg,type);

 Input:
    cfg: a struct, or struct array. Each element in cfg defines 
	 a set of parameters for a simulation. 

    It may contain the following fields:

     *cfg.nphoton:     the total number of photons to be simulated (integer)
     *cfg.prop:        an N by 4 array, each row specifies [mua, mus, g, n] in order.
		       the first row corresponds to medium type 0 which is 
		       typically [0 0 1 1]. The second row is type 1, and so on.
     *cfg.node:        node array for the input tetrahedral mesh, 3 columns: (x,y,z)
     *cfg.elem:        element array for the input tetrahedral mesh, 4 columns
     *cfg.elemprop:    element property index for input tetrahedral mesh
     *cfg.tstart:      starting time of the simulation (in seconds)
     *cfg.tstep:       time-gate width of the simulation (in seconds)
     *cfg.tend:        ending time of the simulation (in second)
     *cfg.srcpos:      a 1 by 3 vector, the position of the source in grid unit
     *cfg.srcdir:      a 1 by 3 vector, specifying the incident vector
     -cfg.facenb:      element face neighbohood list (calculated by faceneighbors())
     -cfg.evol:        element volume (calculated by elemvolume() with iso2mesh)
     -cfg.e0:	       the element ID enclosing the source (calculated by tsearchn())
      cfg.seed:        seed for the random number generator (integer) [0]
      cfg.detpos:      an N by 4 array, each row specifying a detector: [x,y,z,radius]
      cfg.isreflect:   [1]-consider refractive index mismatch, 0-matched index
      cfg.isnormalized:[1]-normalize the output flux to unitary source, 0-no reflection
      cfg.isspecular:  [1]-calculate specular reflection if source is outside
      cfg.basisorder:  [1]-linear basis, 0-piece-wise constant basis
      cfg.outputformat:['ascii'] or 'bin' (in 'double')
      cfg.outputtype:  [X] - output flux, F - fluence, E - energy deposit
      cfg.method:      ray-tracing method, [P]:Plucker, H:Havel (SSE4),
		       B: partial Badouel, S: branchless Badouel (SSE)
      cfg.debuglevel:  debug flag string, a subset of [MCBWDIOXATRPE], no space
      cfg.nout:        [1.0] refractive index for medium type 0 (background)
      cfg.minenergy:   terminate photon when weight less than this level (float) [0.0]
      cfg.roulettesize:[10] size of Russian roulette
      cfg.srctype:     source type, can be ["pencil"],"isotropic" or "cone"
      cfg.srcparam:    1x4 vector for additional source parameter
      cfg.unitinmm:    defines the unit in the input mesh [1.0]

      fields marked with * are required; options in [] are the default values
      fields marked with - are calculated if not given (can be faster if precomputed)

    type: omit or 'omp' for multi-threading version; 'sse' for the SSE4 MMC,
	  the SSE4 version is about 25% faster, but requires newer CPUs.

 Output:
      flux: a struct array, with a length equals to that of cfg.
	    For each element of flux, flux(i).data is a 1D vector with
	    dimensions [size(cfg.node,1) total-time-gates] if cfg.basisorder=1,
	    or [size(cfg.elem,1) total-time-gates] if cfg.basisorder=0. 
	    The content of the array is the normalized flux (or others 
	    depending on cfg.outputtype) at each mesh node and time-gate.
      detphoton: (optional) a struct array, with a length equals to that of cfg.
	    For each element of detphoton, detphoton(i).data is a 2D array with
	    dimensions [size(cfg.prop,1)+1 saved-photon-num]. The first row
	    is the ID(>0) of the detector that captures the photon; the second row
	    is the number of scattering events of the exitting photon; the rest rows
	    are the partial path lengths (in grid unit) traveling in medium 1 up 
	    to the last. If you set cfg.unitinmm, you need to multiply the path-lengths
	    to convert them to mm unit.
      ncfg: (optional), if given, mmclab returns the preprocessed cfg structure,
	    including the calculated subfields (marked by "-"). This can be
	    used as the input to avoid repetitive preprocessing.

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
      % calculate the flux distribution with the given config
      flux=mmclab(cfg);

      cfgs(1)=cfg;
      cfgs(2)=cfg;
      cfgs(1).isreflect=0;
      cfgs(2).isreflect=1;
      cfgs(2).detpos=[30 20 0 1;30 40 0 1;20 30 1 1;40 30 0 1];
      % calculate the flux and partial path lengths for the two configurations
      [fluxs,detps]=mmclab(cfgs);


 This function is part of Mesh-based Monte Carlo (MMC) URL: http://mcx.sf.net/mmc/

 License: GNU General Public License version 3, please read LICENSE.txt for details
</pre>

Because MMCLAB outputs some verbose simulation information, such as 
intermediate timing information, to the stderr descriptor, you will not 
be able to see this information if you launch MATLAB in the GUI mode.
If you want to access to this information, please launch MATLAB with the
following command:

    matlab -nojvm -nodesktop

The debug info will be printed in the console window. Alternatively, you 
can just start the matlab GUI from a console window by typing "matlab",
and the debug information will be printed in the console window instead.


== # Examples ==

We provided several examples to demonstrate the basic usage of MMCLAB,
as well as to perform validations of MMC algorithm using both simple 
homogeneous and heterogeneous domains. These examples are explained below:

==== demo_mcxlab_basic.m ====

In this example, we show the most basic usage of MMCLAB. This include
how to define the input configuration structure, launch MMC simulations
and interpret and plotting the resulting data.

==== demo_validation_homogeneous.m ====

In this example, we validate MMCLAB with a homogeneous medium in a 
cubic domain. This is exactly the example shown in Fig.5 of [Fang2009].

You can also use the alternative optical properties that has a high g
value to observe the similarity between the two scattering/g configurations.

==== demo_validation_heterogeneous.m ====

In this example, we validate the MMCLAB solver with a heterogeneous
domain and the analytical solution of the diffusion model. We also 
demonstrate how to use sub-pixel resolution to refine the representation
of heterogeneities. The domain is consisted of a 6x6x6 cm box with a 
2cm diameter sphere embedded at the center. 

This test is identical to that used for Fig. 3 in [Fang2010].

==== demo_sphere_cube_subpixel.m ====

In this example, we demonstrate how to use sub-pixel resolution 
to represent the problem domain. The domain is consisted of a 
6x6x6 cm box with a 2cm diameter sphere embedded at the center.

==== demo_4layer_head.m ====

In this example, we simulate a 4-layer brain model using MMCLAB.
We will investigate the differences between the solutions with and 
witout boundary reflections (both external and internal) and show
you how to display and analyze the resulting data.


== # How to compile MMCLAB ==

To compile MMCLAB for MATLAB, you need to cd mmc/src directory, and type 

 make mex

from a shell window. You need to make sure your MATLAB is installed and 
the command <tt>mex</tt> is included in your PATH environment variable. Similarly, 
to compile MMCLAB for Octave, you type

 make oct

The command <tt>mkoctfile</tt> must be accessible from your command line
and it is provided in a package named "octave3.x-headers" in Ubuntu (3.x
can be 3.2 or 3.4 etc).


== # Screenshots ==

Screenshot for using MMCLAB in MATLAB:
  http://mcx.sourceforge.net/upload/matlab_mmclab.png

Screenshot for using MMCLAB in GNU Octave:
  http://mcx.sourceforge.net/upload/octave_mmclab.png


== # Reference ==

 [Fang2010] Fang Q, "Mesh-based Monte Carlo method using fast ray-tracing
   in Plucker coordinates," Biomed. Opt. Express 1(1), 165-175 (2010) 

