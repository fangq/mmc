---------------------------------------------------------------------
=                  Mesh-based Monte Carlo (MMC)                     =
        Supporting both OpenCL and Multi-threading with SSE4
---------------------------------------------------------------------


*Author:  Qianqian Fang <q.fang at neu.edu>
*License: GNU General Public License version 3 (GPL v3), see License.txt
*Version: 3.0.0 (v2026, Egg Tart)
*URL:     https://mcx.space/mmc

---------------------------------------------------------------------

Table of Content:

<toc>

---------------------------------------------------------------------

== # What's New ==

MMC v2026 (3.0.0) is a milestone release that brings mesh-based MC to feature parity
with MCX for inverse problems, and makes its ray-tracing verifiable step-by-step.
It is highly recommended to upgrade for all users.

MMC v2026 adds the below key features

* mesh-space adjoint Jacobians (full FEM and nodal approximation) computed on the GPU,
  per-node mua/musp for diffuse optical tomography, and saving both the forward
  solution and the Jacobian in one run
* RF forward and RF/CW adjoint Jacobians ported from MCX, with complex photon weights
  in the mesh-based (BLBadouel) kernel
* MCX-compatible multi-source simulations (`cfg.srcid`, JSON array `Pos`/`Dir`,
  per-slot normalization, `srcid` column in detected photon data)
* new source types: pencilarray, pattern3d, hyperboloid and ring; new OptiX source
  types (disk/gaussian/zgaussian/isotropic/cone/arcsine/line/slit) and a single
  acceleration structure for arbitrary sources in the OptiX-based RT-MMC
  (contributed by Shijie Yan)
* `-F bnii` binary JSON output, `cfg.compileropt`/`cfg.kernelfile` in mmclab/pmmc
* faster CUDA kernel (64 registers, templated RF/multi-source/detphoton paths)
* new functional and step-by-step trajectory tests (`test/testmmc.sh`,
  `test/testtraj.sh`) that run on the CPU, compare MMC with MMC golden data and
  statistically with MCX on matching voxel/mesh domains, and run in the CI
* pmmc is bumped to 0.4.0

this release also fixed a list of bugs, including

* detected photon, trajectory and seed records in `.jdt`/`.jdb` files were tagged
  column-major and were scrambled when loaded by JData readers
* wide-field sources placed on the mesh surface silently deposited the whole beam
  into one element on the GPU; fourierx/fourierx2d source domains were incorrect
* low-mua fluence deposition, detector lookup with extra sources, and OpenCL volume
  buffer allocation (contributed by Lorenzo Pattelli, #129-#132)
* long-standing OpenCL memory leaks (#127), replay with `--optlevel 3` (#122),
  cone/zgaussian launch direction (#121), immc trajectories and partial paths (#116),
  `-F jnii` for detected photons/trajectories (#124), `--maxjumpdebug` (#123),
  unitinmm with `cfg.method='grid'` (#117), faceroi search (#119, #120)

The full changelog is listed below

* 2026-10-21 [ddb5953] [perf] schedule the SSE/OpenMP photon loop dynamically, 1.1-1.3x on 16 threads
* 2026-10-21 [fa723e2] [doc] refresh the v2026 ChangeLog, NEWS and README change lists
* 2026-10-21 [3329d2c] [pmmc] opt-in CUDA backend (-DBUILD_CUDA=ON), fix its device query, test opencl/cuda
* 2026-10-21 [5e5b797] [build] link cudart into the Octave mex built by 'make cudaoct'
* 2026-10-21 [9b7c135] [test] test the opencl, cuda and optix backends compiled into mmc in testmmc.sh
* 2026-10-21 [8b2dfaf] [test] split the trajectory test mesh asymmetrically, regenerate the MMC goldens
* 2026-10-21 [207c6f9] [optix] make the full OptiX debug pipeline opt-in (MMC_OPTIX_DEBUG) instead of !NDEBUG
* 2026-10-21 [1475737] [test] compare trajectories per photon: 5e-2 mm tolerance, tolerate up to 10% forked photons
* 2026-10-21 [ef2c28f] [bug] fix swapped sin/cos in OpenCL MCX_SINCOS with --optlevel 0, rebuild mmc_utils.o when mmc_core.cl changes
* 2026-10-21 [e542dad] [release] set the v2026 (3.0.0) codename to Egg Tart
* 2026-10-21 [64ffdff] [doc] refresh the v2026 ChangeLog, NEWS and README change lists
* 2026-10-21 [de3d1cf] [bug] read a JSON detector object by key regardless of its field count
* 2026-10-21 [eba3e97] [bug] emit valid JSON in the -D S/X/R debug log ("key":[..] instead of "key"=[..])
* 2026-10-21 [91bc95b] [release] bump the remaining v2025.10 keyword strings to v2026
* 2026-10-21 [c1ef008] [doc] v2026 What's New and ChangeLog entry
* 2026-10-21 [e3d94f0] [release] bump version to v2026 (3.0.0), pmmc to 0.4.0, update copyright years
* 2026-10-21 [178598c] [ci] run testmmc.sh, and the trajectory tests on the POCL CPU OpenCL device
* 2026-10-21 [5498290]*[test] add functional, step-by-step trajectory and pmmc tests runnable on the CPU
* 2026-10-21 [2c47938]*[bug] label detected photon, trajectory and seed records in .jdt/.jdb as row-major
* 2026-10-21 [a8f0af4] [bug] mmcsrcdomain: fix fourierx/fourierx2d source domain, treat them as wide-field in mmclab
* 2026-10-21 [757de21]*[bug] mmclab: retessellate wide-field sources by srctype, refuse empty srcelem
* 2026-07-07 [f0b8d82] Avoid volume buffer allocation when OpenCL volume output is disabled
* 2026-07-02 [9f9f1ae] Fix detector lookup offset with extra source records
* 2026-06-29 [d9bc2fe]*Fix low-mua fluence deposition in MMC raytracers
* 2026-06-29 [fc34792] Expose integer isspecular and isextdet in pMMC
* 2026-06-13 [8312b5b] [bug] opencl: harden host buffer handling and fix error-path GPU memory leaks
* 2026-06-13 [ef5bb39]*[bug] fix long-standing memory leakage, due to unused event, fix #127
* 2026-05-21 [5442e3e] [perf] replace NODAL_USE_MUA/MUSP templates with IS_RF/IS_MULTISRC/SAVE_DETPHOTON
* 2026-05-21 [f8c8cb7] [perf] cuda: __launch_bounds__(64,16) on mmc_main_loop drops 122->64 regs
* 2026-05-20 [9dc5ec0] [clean] mmc: revert per-slot energy accumulator to single-bin + broadcast
* 2026-05-20 [12652da] [bug] mesh adjoint multi-source: per-slot initial tet + per-slot energy accounting, w/ Claude code
* 2026-05-14 [35cfefe] [bug] pmmc: build detector-as-adjoint srcdata BEFORE mesh_validate, w/ Claude code
* 2026-05-14 [47e5406] [feat] pmmc: accept Nsrc-row srcpos / srcdir for multi-source forward, w/ Claude code
* 2026-05-13 [6a90384] [clean] suppress OpenCL build log on successful kernel build, w/ Claude code
* 2026-05-13 [f769204] [fix] silence OpenCL JIT warnings on srcid multi-source guard, w/ Claude code
* 2026-05-13 [30498af] [feat] mmclab honors USE_MCXCL=N device override, ported from mcxlabcl, w/ Claude code
* 2026-05-13 [23cde0a] [feat] mmclab/pmmc adjoint detector append preserves multi-source srcdata, w/ Claude code
* 2026-05-13 [4b26968] [format] apply astyle + mh_style fixes from make pretty / mh_style --fix .
* 2026-05-13 [4e6550e] [fix] guard sgtitle() in mmclab demos for Octave < 8.0 compat, w/ Claude code
* 2026-05-13 [8868080]*[feat] RF complex weight tracking in BLBadouel mesh kernel + complex flux.data in mmclab mesh mode, w/ Claude code
* 2026-05-13 [461b1ed] [feat] mmclab: multi-source / mesh adjoint demos + 2D srcpos parser + RF imag for mesh basisorder=1, w/ Claude code
* 2026-05-13 [abeb635] [feat] per-slot normalizer for extra-source slots scales by srcdata weight ratio, w/ Claude code
* 2026-05-13 [a90ef05]*[feat] mcx-parity multi-source: JSON array Pos/Dir, srcid>0 single-slot, detphoton srcid column, w/ Claude code
* 2026-05-13 [c39d979] [bug] fix per-node mua/musp parsers + missing OpenCL JIT macros, w/ Claude code
* 2026-05-12 [5433129]*[feat] per-node mua/musp via CUDA templates + OpenCL macros for DOT recon, w/ Claude code
* 2026-05-12 [0d8fb36] [bug] fix mesh-adjoint flux.data layout + nsrcslots gate, w/ Claude code
* 2026-05-12 [ed4733e]*[feat] mesh-space adjoint Jacobian (full FEM + nodal approx) via GPU kernels, w/ Claude code
* 2026-05-12 [959562b] [bug] fix maxgate>1 adjoint buffer collision: framelen should be datalen, w/ Claude code
* 2026-05-12 [25b5982]*[feat] add cfg.srcid and mesh-mode multi-source slot infrastructure, w/ Claude code
* 2026-05-10 [76f945d]*[feat] support -F bnii to save binary json files, use -Z base64 for raw, faster IO
* 2026-05-10 [2a7c1f8] [bug] enable -F jnii flag for saving detphoton and traj data, fix #124
* 2026-05-10 [31d6a06] [bug] add the missing --maxjumpdebug flag to cli, fix #123
* 2026-05-08 [82994b1] [ci] add the missing detparam1/detparam2 input interface in JSON
* 2026-04-04 [47b8472]*[feat] save and export both forward model and the Jacobian outputs
* 2026-04-04 [2434e3a]*[feat] port RF forward and adjoint Jacobian from mcx to mmc
* 2026-03-27 [a1465d4] [bug] fix additional seed buffer bugs to enable replay, fix #122, w/ Claude code
* 2026-03-27 [90081a3] [bug] fix replay failure when setting optlevel>=3, fix #122, bump pmmc 0.3.9
* 2026-03-27 [c2ea9e1] [bug] fix replay failure when setting optlevel>=3, fix #122, bump pmmc 0.3.9
* 2026-03-22 [db6a0ed]*[feat] support pencilarray, pattern3d, hyperboloid and ring src, must retesselate
* 2026-03-22 [fd74cbc] [feat] support cfg.compileropt and cfg.kernelfile in mmclab and pmmc
* 2026-03-22 [f4052eb] [bug] fix mmc cone beam and zgaussian launch dir rotation, fix #121
* 2026-03-06 [4cbf7b1] [bug] force input elem to mmcaddsrc and mmcadddet to have label
* 2026-01-11 [1799480] [bug] stop testing faceroi in the 1st neighbor, close #119, close #120
* 2025-12-31 [47617fd] [ci] remove octave image from mcxcl package
* 2025-12-14 [85a7bfa] [bug] debug immc ray-tracing, fix incorrect partial path output, #116
* 2025-12-14 [f03d523] [debug] change debug ouput to json, fix cfg.debugphoton bug
* 2025-12-10 [95dc994] [ci] remove retired macos-13 runner, use macos-15-intel
* 2025-12-10 [85613f9] [bug] fix missing trajectories in the immc mode, #116
* 2025-12-10 [2a58c43] [bug] fix disabled unitinmm setting when cfg.method='grid', fix #117
* 2025-12-01 [0eedf2b]*[optix] add source types: disk/gaussian/zgaussian/isotropic/cone/arcsine/line/slot
* 2025-12-01 [3fe4024] [optix] fix some memory leakage
* 2025-11-30 [639e810] [optix] add back normalization, at the moment, energytot is set to nphoton
* 2025-11-30 [f048bd9] [optix] fix regression due to float3/FLOAT3 transition, patched by Shijie
* 2025-11-30 [263d23f] [optix] fix example script path
* 2025-11-30 [5b639c6] [format] automatically reformat all new matlab .m scripts
* 2025-11-30 [72015ad] [optix] fix ifdef in mmclab
* 2025-11-29 [3dd250a] [optix] reformat all optix units, rename all units start with mmc_
* 2025-11-28 [5402644] [ci] fix g++-12 missing on macos-14
* 2025-11-28 [deb927f] [ci] fix g++-9 missing on ubuntu 22.04
* 2025-11-28 [7be94c0] [bug] fix immc unable to save ROI partial-path, fix #116, pmmc 0.3.8
* 2025-11-17 [ad54567] remove redundant debugging code
* 2025-11-17 [9ff1501] fix tissue labels for wide-field simulations
* 2025-11-17 [5e52cb2] add planar source in RT-DMMCLAB
* 2025-11-16 [fada4bc] search init. medium ID at runtime, add medium ID for terminating photons
* 2025-11-13 [5e5cb81] [pmmc] bump pmmc version to 0.3.7 after fix #108
* 2025-11-13 [e55aa08] [bug] allow passing srcdir(4) as NaN and -Inf, accept --json, fix #108
* 2025-11-10 [f65b5a7] [doc] explain widefield detector settings in mmclab help info
* 2025-10-25 [4ebbcbb] [mem] optimize struct align param to minimize padding
* 2025-10-25 [563aaeb] [std] make unit compatible with some compilers
* 2023-10-27 [d47862c] print kernel speed
* 2023-10-27 [4cd86ca] improve variable definition
* 2023-06-09 [acefed0] add slen to optix trace call
* 2022-09-28 [2571307] disable anyhit program to improve speed
* 2022-09-22 [1f4608f] free the host seed buffer after use
* 2022-09-22 [47126c4] update example to output normalized fluence
* 2022-09-21 [c998f72] add double buffer for OpenCL and CUDA MMC to reduce numerical loss
* 2022-09-17 [9891182] add output normalization for OptiX
* 2022-09-17 [6a40a5b] add fluence output in OptiX
* 2022-09-15 [fc45dd9] add reflection/refraction for mismatched boundaries
* 2022-09-14 [fb50bf6] add face normal vector to SBT
* 2022-09-14 [1dd44e3] reduce miss cases by setting tmin to 0.0f
* 2022-09-10 [5ba4606] add debug messages for context creation
* 2022-09-08 [14fa090] no need to check volume mesh for triangle-based MC
* 2022-09-07 [10dc9ed] fix lower 16-bit mask, thanks to Douglas!
* 2022-09-07 [5c1ee41] fix an edge case for matched media bounaries
* 2022-09-07 [834669b] remove redundant compilation files
* 2022-09-07 [ec74b75] automate compilation for optix dependencies
* 2022-09-07 [9f48be7] allow GNU make to track changes in header files
* 2022-09-06 [2a6bfa8] resolve the mex compiler complaint
* 2022-09-06 [b6717c6] save output only when moving to a new element
* 2022-09-06 [6168b64] Reorganize code for better readability
* 2022-09-06 [208027a] tweak tmin threshold to work with colin27 example
* 2022-09-06 [4eb13c6] add function to get surface mesh from volume mesh, thanks to Qianqian!
* 2022-09-05 [d477bc5] add command line options and fix bugs
* 2022-09-05 [86d1686] add triangular surface-based MC



------------------------------------------------------------------------------- 

== # Introduction ==

Mesh-based Monte Carlo (MMC) is a 3D Monte Carlo (MC) simulation software 
for photon transport in complex turbid media. MMC combines the strengths
of the MC-based technique and the finite-element (FE) method: on the 
one hand, it can handle general media, including low-scattering ones, 
as in the MC method; on the other hand, it can use an FE-like tetrahedral 
mesh to represent curved boundaries and complex structures, making it
even more accurate, flexible, and memory efficient. MMC uses the
state-of-the-art ray-tracing techniques to simulate photon propagation in 
a mesh space. It has been extensively optimized for excellent computational
efficiency and portability. MMC currently supports multi-threaded 
parallel computing via OpenMP, Single Instruction Multiple Data (SIMD) 
parallelism via SSE and, starting from v2019.10, OpenCL to support a wide
range of CPUs/GPUs from nearly all vendors.

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

The most relevant publication describing this work is the GPU-accelerated
MMC paper:

Qianqian Fang and Shijie Yan, "GPU-accelerated mesh-based Monte Carlo 
photon transport simulations," J. of Biomedical Optics, in press, 2019.
Preprint URL: https://www.biorxiv.org/content/10.1101/815977v1

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

== # Download and Compile MMC ==

The latest release of MMC can be downloaded from the following URL:

  https://mcx.space/#mmc

The development branch (not fully tested) of the code can be accessed 
using Git. However this is not encouraged unless you are
a developer. To check out the Git source code, you should use the following 
command:

  git clone https://github.com/fangq/mmc.git mmc

To compile the software, you need to install GNU gcc compiler toolchain
on your system. For Debian/Ubuntu based GNU/Linux systems, you can type

  sudo apt-get install build-essential

and for Fedora/Redhat based GNU/Linux systems, you can type

  sudo dnf install make automake gcc gcc-c++
 
To compile the binary with multi-threaded computing via OpenMP, 
your `gcc` version should be at least 4.0. To compile the binary 
supporting SSE4 instructions, gcc version should be at least 
4.3.4. For windows users, you should install MSYS2 or Cygwin64 [3]. During the 
installation, please select `mingw64-x86_64-gcc` and `make` packages.
For Mac OS users, you need to install gcc from 
Homebrew or MacPorts and use the instructions below to compile the 
MMC source code.

=== # Building MMC using CMake ===

One can choose one of the two methods to build mmc binaries. The first
approach is to use CMake. CMake is a portable system creating compilation
and linking commands automatically adapted to your operating system and 
installed compiler. It can run on Linux, MacOS and Windows.

To use CMake, you will have to first run `sudo apt-get install cmake`
or `sudo dnf install cmake` to install cmake first. To build MMC binaries,
you should first navigate to the `mmc/src` folder, and run

  mkdir -p build
  cd build
  cmake .. && make

if cmake complains that any required library is missing, you will need
to install those dependencies, removing all files inside the build folder,
and run the cmake command above again.

The above command builds the `mmc` executable inside `mmc/bin` folder. If
your system has MATLAB installed, the above command also builds mmclab mex
file as `mmc/mmclab/mmc.mex*` where the mex suffix depends on your OS.

If you want to build the "Trinity" version of mmc to support CUDA on NVIDIA
GPUs, you will have to first install CUDA toolkit, and replace the above cmake
command by 

  cmake .. -DBUILD_CUDA=on && make

the executable will be build as `mmc/bin/mmciii` and the mex file is `mmc/mmclab/mmciii.mex*`.


=== # Building MMC using GNU Make ===

To compile the program, you should first navigate into the mmc/src folder,
and type

  make

this will create a fully optimized OpenCL based mmc executable, 
located under the mmc/src/bin/ folder.

Other compilation options include

  make ssemath  # this is the same as make, building mmc binary with SSE4+OpenMP+OpenCL
  make cuda     # this compiles the "Trinity" version of mmc, supports SSE4+OpenMP+OpenCL+CUDA
  make omp      # this compiles a multi-threaded binary using OpenMP
  make release  # create a single-threaded optimized binary
  make prof     # this makes a binary to produce profiling info for gprof
  make sse      # this uses SSE4 for all vector operations (dot, cross), implies omp

if you wish to build the mmc mex file to be used in matlab, you should run

  make mex      # this produces mmc.mex* under mmc/mmclab/ folder
  make cudamex  # this produces a "Trinity" version of mmc.mex* that supports SSE+OpenCL+CUDA

similarly, if you wish to build the mex file for GNU Octave, you should run

  make oct      # this produces mmc.mex* under mmc/mmclab/ folder
  make cudaoct  # this produces a "Trinity" version of mmc.mex* that supports SSE+OpenCL+CUDA

If you append `-f makefile_sfmt` at the end of any of the above 
make commands, you will get an executable named `mmc_sfmt`, which uses a 
fast MT19937 random-number-generator (RNG) instead of the default GLIBC 
48bit RNG. If your CPU supports SSE4, the fastest binary can be obtained
by running the following command:

  make ssemath -f makefile_sfmt

You should be able to compile the code with an Intel C++ compiler,
an AMD C compiler or LLVM compiler without any difficulty. To use other
compilers, you simply append `CC=compiler_exe` to the above make 
commands. If you see any error messages, please google and fix 
your compiler settings or install the missing libraries.

A special note for Mac OS users: you can you use both gcc (installed by MacPorts
or brew) or the default clang gcc provided by the default Xcode compiler to build
mmc. MMC requires OpenMP for multi-threading based parallel computing. If one uses
the clang compiler, one must first install `libomp` package in order to compile mmc.

  brew install libomp
  brew link --force libomp

One can switch to other compilers by setting the `CC`, `CXX` and `AR` environment
variables, for example

  make CC=gcc-mp-10 CXX=g++-mp-10 AR=g++-mp-10

After compilation, you may add the path to the "mmc" binary (typically,
mmc/src/bin) to your search path. To do so, you should modify your 
$PATH environment variable. Detailed instructions can be found at [5].

You can also compile MMC using Intel's C++ compiler - icc. To do this, you run

  make CC=icc

you must enable `icc` related environment variables by source the `compilervars.sh`
file. The speed of icc-generated mmc binary is generally faster for CPU/SSE based
MMC simulation than those compiled by `gcc`.

-------------------------------------------------------------------------------

== # Running Simulations ==


=== Preparation ===

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


=== Command line options ===

The full command line options of MMC include the following:
<pre>
###############################################################################
#                     Mesh-based Monte Carlo (MMC) - OpenCL                   #
#          Copyright (c) 2010-2025 Qianqian Fang <q.fang at neu.edu>          #
#              https://mcx.space/#mmc  &  https://neurojson.io/               #
#                                                                             #
#Computational Optics & Translational Imaging (COTI) Lab  [http://fanglab.org]#
#   Department of Bioengineering, Northeastern University, Boston, MA, USA    #
###############################################################################
#    The MCX Project is funded by the NIH/NIGMS under grant R01-GM114365      #
###############################################################################
#  Open-source codes and reusable scientific data are essential for research, #
# MCX proudly developed human-readable JSON-based data formats for easy reuse.#
#                                                                             #
#Please visit our free scientific data sharing portal at https://neurojson.io/#
# and consider sharing your public datasets in standardized JSON/JData format #
###############################################################################
$Rev::      $v2026$Date::                       $ by $Author::             $
###############################################################################

usage: mmc <param1> <param2> ...
where possible parameters include (the first item in [] is the default value)

== Required option ==
 -f config     (--input)       read an input file in .json or inp format

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
 -M [G|SG]    (--method)      choose ray-tracing algorithm (only use 1 letter)
                               P - Plucker-coordinate ray-tracing algorithm
                               H - Havel's SSE4 ray-tracing algorithm
                               B - partial Badouel's method (used by TIM-OS)
                               S - branch-less Badouel's method with SSE
                               G - dual-grid MMC (DMMC) with voxel data output
 -e [1e-6|float](--minenergy)  minimum energy level to trigger Russian roulette
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh
 -k [1|0]      (--voidtime)    when src is outside, 1 enables timer inside void

== GPU options ==
 -A [0|int]    (--autopilot)   auto thread config:1 enable;0 disable
 -c [opencl,sse,cuda](--compute) select compute backend (default to opencl)
                               can also use 0: sse, 1: opencl, 2: cuda
 -G [0|int]    (--gpu)         specify which GPU to use, list GPU by -L; 0 auto
      or                       if set to -1, CPU-based SSE mmc will be used
 -G '1101'     (--gpu)         using multiple devices (1 enable, 0 disable)
 -W '50,30,20' (--workload)    workload for active devices; normalized by sum
 --atomic [1|0]                1 use atomic operations, 0 use non-atomic ones

== Output options ==
 -s sessionid  (--session)     a string used to tag all output file names
 -O [X|XFEJLP] (--outputtype)  X - output flux, F - fluence, E - energy density
                               J - Jacobian, L - weighted path length, P -
                               weighted scattering count (J,L,P: replay mode)
 -d [0|1]      (--savedet)     1 to save photon info at detectors,0 not to save
 -H [1000000] (--maxdetphoton) max number of detected photons
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
 -F [bin|...] (--outputformat) 'ascii', 'bin' (in 'double'), 'mc2' (double) 
                               'hdr' (Analyze) or 'nii' (nifti, double)
                               mc2 - MCX mc2 format (binary 64bit float)
                               jnii - JNIfTI format (https://neurojson.org)
                               bnii - Binary JNIfTI (https://neurojson.org)
                               nii - NIfTI format
                               hdr - Analyze 7.5 hdr/img format
    the bnii/jnii formats support compression (-Z) and generate small files
    load jnii (JSON) and bnii (UBJSON) files using below lightweight libs:
      MATLAB/Octave: JNIfTI toolbox   https://github.com/NeuroJSON/jnifti, 
      MATLAB/Octave: JSONLab toolbox  https://github.com/fangq/jsonlab, 
      Python:        PyJData:         https://pypi.org/project/jdata
      JavaScript:    JSData:          https://github.com/NeuroJSON/jsdata
 -Z [zlib|...] (--zip)      set compression method if -F jnii or --dumpjson
                            is used (when saving data to JSON/JNIfTI format)
                            0 zlib: zip format (moderate compression,fast) 
                            1 gzip: gzip format (compatible with *.gz)
                            2 base64: base64 encoding with no compression
                            3 lzip: lzip format (high compression,very slow)
                            4 lzma: lzma format (high compression,very slow)
                            5 lz4: LZ4 format (low compression,extrem. fast)
                            6 lz4hc: LZ4HC format (moderate compression,fast)
 --dumpjson [-,2,'file.json'] export all settings, including volume data using
                          JSON/JData (https://neurojson.org) format for 
                          easy sharing; can be reused using -f
                          if followed by nothing or '-', mmc will print
                          the JSON to the console; write to a file if file
                          name is specified; by default, prints settings
                          after pre-processing; '--dumpjson 2' prints 
                          raw inputs before pre-processing

== User IO options ==
 -h            (--help)        print this message
 -v            (--version)     print MMC version information
 -l            (--log)         print messages to a log file instead
 -i            (--interactive) interactive mode

== Debug options ==
 -D [0|int]    (--debug)       print debug information (you can use an integer
  or                           or a string by combining the following flags)
 -D [''|SCBWDIOXATRPEM]        1 S  photon movement info
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
                            8192 M  return photon trajectories
      combine multiple items by using a string, or add selected numbers together
 --debugphoton [-1|int]        to print the debug info specified by -D only for
                               a single photon, followed by its index (start 0)

== Additional options ==
 --momentum     [0|1]          1 to save photon momentum transfer,0 not to save
 --gridsize     [1|float]      if -M G is used, this sets the grid size in mm
 --maxjumpdebug [10000000|int] when trajectory is requested (i.e. -D S),
                               use this parameter to set the maximum positions
                               stored (default: 1e7)

== Example ==
       mmc -n 1000000 -f input.json -s test -b 0 -D TP -G -1
</pre>


=== Input files ===

It is highly recommended to use the JSON-formatted input file described in the following
section. The legacy input file format `.inp` is depreciated and may be removed in
future releases.

The simplest example can be found under the "example/onecube" 
folder. Please run "createmesh.m" first from Matlab/Octave to 
create all the mesh files, which include
<pre>
elem_onecube.dat    -- tetrahedral element file
facenb_onecube.dat  -- element neighbors of each face
node_onecube.dat    -- node coordinates
prop_onecube.dat    -- optical properties of each element type
velem_onecube.dat   -- volume of each element
</pre>

The input file of the example is named "onecube.inp", where we
specify most of the simulation parameters. The input file follows
a similar format as in MCX, which looks like the following

<pre>
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
</pre>

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


=== JSON-formatted input files ===

Starting from version 0.9, MMC accepts a JSON-formatted input file in
addition to the conventional tMCimg-like input format. JSON 
(JavaScript Object Notation) is a portable, human-readable and 
"fat-free" text format to represent complex and hierarchical data.
Using the JSON format makes a input file self-explanatory, extensible
and easy-to-interface with other applications (like MATLAB).

A sample JSON input file can be found under the examples/onecube
folder. The same file, onecube.json, is also shown below:
<pre>
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
</pre>

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
<pre>
{
    "Domain.MeshID": "onecube",
    "Session.ID": "onecube",
    ...
}
</pre>

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


=== Photon debugging information using -D flag ===

the output format for -D M (photon moving) is below:
<pre>
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
</pre>

-------------------------------------------------------------------------------

=== Plotting the Results ===

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

== # Known issues and TODOs ==

* MMC only supports linear tetrahedral elements at this point. Quadratic \
 elements will be added later
* Currently, this code only supports element-based optical properties; \
 nodal-based optical properties (for continuously varying media) will be \
 added in a future release

-------------------------------------------------------------------------------

== # Getting Involved ==

MMC is an open-source software. It is released under the terms of GNU 
General Public License version 3 (GPLv3). That means not only everyone 
can download and use MMC for any purposes, but also you can modify the 
code and share the improved software with others (as long as the derived 
work is also licensed under the GPLv3 license).

If you already made a change to the source code to fix a bug you encountered 
in your research, we are appreciated if you can share your changes (as 
`git diff` outputs) with the developers. We will patch the code as soon 
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

 https://mcx.space/#about

-------------------------------------------------------------------------------

== # Acknowledgement ==

MMC uses the following open-source libraries:

=== ZMat data compression unit ===

* Files: src/zmat/*
* Copyright: 2019-2020 Qianqian Fang
* URL: https://github.com/fangq/zmat
* License: GPL version 3 or later, https://github.com/fangq/zmat/blob/master/LICENSE.txt

=== LZ4 data compression library ===

* Files: src/zmat/lz4/*
* Copyright: 2011-2020, Yann Collet
* URL: https://github.com/lz4/lz4
* License: BSD-2-clause, https://github.com/lz4/lz4/blob/dev/lib/LICENSE

=== LZMA/Easylzma data compression library ===

* Files: src/zmat/easylzma/*
* Copyright: 2009, Lloyd Hilaiel, 2008, Igor Pavlov
* License: public-domain
* Comment: \
 All the cruft you find here is public domain.  You don't have to \
 credit anyone to use this code, but my personal request is that you mention \
 Igor Pavlov for his hard, high quality work.

=== Miniz compression library ===

- Files: src/zmat/miniz/*
- Copyright 2013-2014 RAD Game Tools and Valve Software
- Copyright 2010-2014 Rich Geldreich and Tenacious Software LLC
- URL: https://github.com/richgel999/miniz
- License: MIT-license, https://github.com/richgel999/miniz/blob/master/LICENSE

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

* Redistributions of source code must retain the above copyright \
 notice, this list of conditions and the following disclaimer. 
* Redistributions in binary form must reproduce the above \
 copyright notice, this list of conditions and the following \
 disclaimer in the documentation and/or other materials provided \
 with the distribution.
* Neither the name of the Hiroshima University nor the names of \
 its contributors may be used to endorse or promote products \
 derived from this software without specific prior written \
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

== # Reference ===

[1] http://iso2mesh.sf.net  -- an image-based surface/volumetric mesh generator
[2] https://mcx.space       -- Monte Carlo eXtreme: a GPU-accelerated MC code
[3] https://cygwin.com/setup-x86_64.exe
[4] http://developer.apple.com/mac/library/releasenotes/DeveloperTools/RN-llvm-gcc/index.html
[5] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/AddPath
[6] https://mcx.space/cgi-bin/index.cgi?MMC/Doc/FAQ#How_do_I_interpret_MMC_s_output_data
[7] http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?fun/qmeshcut
