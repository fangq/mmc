---------------------------------------------------------------------
=                    Mesh-based Monte Carlo (MMC)                   =
---------------------------------------------------------------------

*Author: Qianqian Fang (q.fang at neu.edu)
*License: GNU General Public License version 3 (GPLv3)
*Version: 3.0.0 (v2026, Egg Tart)
*Website: https://mcx.space/#mmc
*Download: https://mcx.space/wiki/?Get
*Discussion: https://groups.google.com/g/mmc-users (MMC) and https://groups.google.com/g/mcx-users (all MCX software)
---------------------------------------------------------------------

Table of Content:

<toc>

---------------------------------------------------------------------

MMC is a Monte Carlo simulator of light transport in 3D turbid media, such as
biological tissues, represented by '''tetrahedral meshes'''. Meshes follow curved
boundaries and thin layers far more accurately than voxels, and use fewer elements
where the anatomy is simple. MMC runs on multi-core CPUs and on GPUs through four
interchangeable backends - SSE/OpenMP, OpenCL, CUDA and OptiX - and shares its input
format, MATLAB/Python interfaces and output formats with its voxel-based sibling,
[https://github.com/fangq/mcx MCX].

MMC can be used in four ways:

<pre>
`mmc` command line
    What it is: a standalone executable reading a JSON input file
    Best for: batch jobs, clusters, scripting, reproducible inputs
MMCLAB (MATLAB/Octave)
    What it is: `mmclab(cfg)`, which also prepares the mesh for you
    Best for: meshing and simulating in one script, image reconstruction
PMMC (Python)
    What it is: `pmmc.run(cfg)` / `pmmc.mmclab(cfg)`, installed with `pip`
    Best for: Python/NumPy workflows with pyiso2mesh (https://pypi.org/project/iso2mesh/)
MCX Studio (GUI)
    What it is: a graphical editor and viewer that runs mcx, mcxcl and mmc
    Best for: inspecting meshes, sources and results
</pre>

== # What's New ==

MMC v2026 (3.0.0) is a milestone release that brings mesh-based MC to feature parity
with MCX for inverse problems, and makes its ray-tracing verifiable step-by-step.
It is highly recommended to upgrade for all users.

MMC v2026 adds the below key features

* directly building adjoint Jacobians via `cfg.outputtype='adjoint'`: mesh-space
  Jacobians (full FEM and nodal approximation) computed on the GPU, per-node mua/musp
  for diffuse optical tomography, and saving both the forward solution and the
  Jacobian in one run
* frequency-domain (FD/RF) forward simulations and adjoint-based FD-Jacobians
  (complex-valued), ported from MCX, with complex photon weights in the mesh-based
  (BLBadouel) kernel
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
* faster simulations: the SSE/OpenMP photon loop is scheduled dynamically (1.1-1.3x
  on 16 threads), the OpenCL kernel is built with `-cl-unsafe-math-optimizations`
  (1.14-1.16x) and photons are dispatched dynamically in the OpenCL/CUDA kernels
  (1.02-1.06x); because of the dynamic scheduling, multi-threaded runs with the same
  seed agree statistically but are not bit-reproducible, use a single thread
  (`OMP_NUM_THREADS=1`) on the CPU to reproduce a run exactly
* the OptiX-based RT-MMC now supports all source types, multiple sources, detected
  photons, trajectories (`-D M`), adjoint Jacobians and surface-only meshes
* implicit MMC (node/edge/face ROIs and `Mesh.MeshROIPrim`) now also runs on the
  OpenCL and CUDA backends, ported from the GPU iMMC of Aiden Lewis: each element
  holds a host-built list of all overlapping capsules and spheres, traced as a union,
  and face slabs apply Fresnel reflection on single-sided slabs; the GPU results agree
  with the CPU within 0.1% and pass the analytic iMMC tests
* mmc-optix packages for Linux and Windows: mmc with the CUDA (`-c cuda`) and OptiX
  (`-c optix`) backends next to the CPU and OpenCL ones, for NVIDIA GPUs with driver
  R515 or newer; mmc now also builds with MSVC
* the adjoint Jacobians of mus/musp apply the same chain rule as MCX
* `pmmc.mmclab()` applies the `mmclab.m` pre/post-processing around `pmmc.run()`, and
  all mmclab example scripts are ported to Python under `pmmc/example`
* unit tests and coverage: host-side unit tests of the parsers, mesh code, ray tracers
  and the CUDA/OptiX backends (`make unittest`, 513 checks) raise the line coverage of
  the host code to about 95%; 68 built-in example tests check every run against golden
  values; the bugs these tests found were fixed
* fuzzing: rounds of JSON-aware mutation fuzzing of the input files and of the
  command-line options, run on AddressSanitizer builds of mcx, mcxcl and mmc (for mmc,
  `test/fuzz/mmcfuzz.py`, about 8,000 cases on the CPU, OpenCL, CUDA and OptiX backends); every
  crash, hang and memory error they found was fixed
* hardening: JSON arrays are length-checked before use; JData arrays must carry valid
  annotations and decode to their declared size; an embedded mesh is checked (node
  indices, element labels, finite coordinates, InitElem); optical properties, time
  gates, SrcNum and pattern sizes are validated; string options are bounds-checked;
  invalid input now stops with an error message instead of crashing or hanging
* fixes found by fuzzing: photons reaching the time limit on the CPU path are stopped as
  on the GPU instead of being retraced with an undefined element; the last step before
  the time limit can no longer be negative; lost photons no longer read outside the
  mesh; pattern sizes are bounded before their use as indices
* frequency-domain (RF) forward simulations now also run on the CPU (`-c sse`), and RF
  output is fixed for element- and node-based runs on all backends (a 19-degree phase
  error, and a missing imaginary part on the GPU); the standalone binary saves the
  imaginary part to `<session>_imag.<ext>`; CW, RF and time-domain solutions of every
  backend were checked against the diffusion solution of a semi-infinite medium (within
  2% in amplitude and 0.3 degrees in phase)
* `Session.OutputType` in JSON inputs accepts the names mmclab and pmmc use; only the
  first letter used to be read
* the released Linux binaries and mex files are built in a manylinux2014 container and
  run on glibc 2.17 or newer; Octave mex files for Windows and macOS are built against
  Octave 10
* pmmc is bumped to 0.3.12

this release also fixed a list of bugs, including

* the long-standing OpenCL GPU memory leak, first reported in 2020, is finally fixed
  (#127): on NVIDIA's OpenCL runtime, every simulation lost about 300-400 MB of GPU
  memory, so repeated mmclab/pmmc calls in one MATLAB/Octave/Python session eventually
  ran out of GPU memory; an unreleased kernel-launch event kept a reference to the
  OpenCL context, so the GPU buffers were never freed. AMD, Intel and CPU (POCL)
  devices were not affected
* detected photon, trajectory and seed records in `.jdt`/`.jdb` files were tagged
  column-major and were scrambled when loaded by JData readers
* wide-field sources placed on the mesh surface silently deposited the whole beam
  into one element on the GPU; fourierx/fourierx2d source domains were incorrect
* low-mua fluence deposition, detector lookup with extra sources, and OpenCL volume
  buffer allocation (contributed by Lorenzo Pattelli, #129-#132)
* replay with `--optlevel 3` (#122), cone/zgaussian launch direction (#121), immc
  trajectories and partial paths (#116),
  `-F jnii` for detected photons/trajectories (#124), `--maxjumpdebug` (#123),
  unitinmm with `cfg.method='grid'` (#117), faceroi search (#119, #120)
* no fluence is deposited in non-tissue (label<=0) elements, and the nodal output
  uses the tissue-only nodal volume consistently on all backends
* the scattering angle and the albedo weight in implicit MMC (edge/node/face ROIs)
  now use the optical properties of the ROI instead of the host element
* the host-generated RNG seeds were not uniformly random (bit 31 was always 0 with
  the 31-bit rand() of Linux/macOS, only 15 bits per word with MSVC); seeds are now
  uniform on all backends, which changes the random streams for a given seed
* OptiX: no fluence in the ambient medium, rays that exit the mesh or miss all
  triangles are terminated, and OptiX builds without the SDK's `sutil` headers
* large meshes are no longer rejected by the CUDA backend (off-by-one MAX_PROP check)
* pmmc macOS wheels are tagged with their actual architecture instead of
  `universal2`, and the Intel macOS module links the OpenMP runtime

The most important changes since v2025.10 are summarized below; `ChangeLog.txt`
lists all of them, with their commit hashes.

'''Inverse problems and frequency domain'''
* [ed4733e] mesh-space adjoint Jacobians (full FEM and nodal approximation) computed on the GPU;
  [14bf357] mus/mus' Jacobians use the same chain rule as MCX
* [5433129] per-node mua/mus' (`cfg.nodemua`, `cfg.nodemusp`) for diffuse optical tomography
* [2434e3a] RF forward and RF/CW adjoint Jacobians ported from MCX; [8868080] complex photon
  weights in the mesh kernel; [c87499f] RF forward on the CPU, and RF output fixed for element-
  and node-based runs
* [47b8472] the forward solution and the Jacobians are saved together

'''Sources, multi-source runs and outputs'''
* [25b5982, a90ef05] MCX-compatible multi-source simulations (`cfg.srcid`, JSON arrays of
  `Pos`/`Dir`, a srcid column in detected photon data)
* [db6a0ed] pencilarray, pattern3d, hyperboloid and ring sources
* [76f945d] `-F bnii` binary JSON output; [621f207] parallel zlib/gzip compression of JData data
* [fd74cbc] `cfg.compileropt` and `cfg.kernelfile` in mmclab and pmmc

'''Implicit MMC'''
* [376ea5d] iMMC traces the union of all capped capsules and spheres overlapping an element;
  robust slab tests; [a2b303c] iMMC on the OpenCL and CUDA backends (ported from Aiden Lewis'
  GPU iMMC); [a4a3f6c] `Mesh.MeshROIPrim` capsules and spheres given by coordinates
* [5969413, f6cf915] faster GPU kernels: element normals in global memory, one loop of ray-tet steps

'''OptiX (RT-MMC) and CUDA'''
* [0eedf2b, edaaa02] all source types in OptiX MMC (contributed by Shijie Yan); [2434293] multiple
  sources; [7ee0a9d] detected photons; [6ac3389] trajectories; [6b92686] adjoint Jacobians;
  [3671f6d] surface-only meshes; [0ec2400] unsupported settings are rejected with a message
* [f8c8cb7] CUDA kernel down from 122 to 64 registers
* [24aebad] mmc builds with MSVC, enabling CUDA/OptiX on Windows; [bfff982] mmc-optix packages
  for Linux and Windows

'''Performance'''
* [ddb5953] dynamic OpenMP scheduling on the CPU, 1.1-1.3x on 16 threads
* [db4d829] OpenCL kernel built with `-cl-unsafe-math-optimizations`, 1.14-1.16x
* [e803195] dynamic photon dispatch in the OpenCL/CUDA kernels, 1.02-1.06x

'''Bug fixes'''
* [2c47938] detected photon, trajectory and seed records in `.jdt`/`.jdb` were scrambled
  when loaded by JData readers
* [757de21] wide-field sources on the mesh surface deposited the whole beam into one element on
  the GPU; fourierx/fourierx2d source domains fixed
* [a5e9c1f] the GPU lost up to 40% of the diffuse reflectance
* [ef5bb39] long-standing OpenCL memory leaks (#127)
* [f6091aa] host RNG seeds were not uniformly random; now uniform on all backends
* [b64aa2a] no fluence in non-tissue (label<=0) elements; consistent nodal volumes
* [882054e] implicit MMC uses the ROI's optical properties for scattering and albedo
* [66f4c90] `cfg.face` is ignored when `cfg.elem` is given (fixing redbird's MMC Jacobians)
* [816f671] `Session.OutputType` names read in full, not just their first letter

'''MMCLAB and PMMC'''
* [5ecd5ac] `pmmc.mmclab()`: the mmclab.m pre/post-processing around `pmmc.run()`; [d820fb5] all
  mmclab examples ported to Python under `pmmc/example`
* [a6d4f42] the mmclab help rewritten to document every cfg field; [0feb5bb] `cfg.issave2pt`,
  and an error for backends the mex was not built with

'''Testing and robustness'''
* [5498290] functional, step-by-step trajectory and pmmc tests that run on the CPU in CI
* [b07e548, 94de589] host unit tests, including the CUDA/OptiX backends; host coverage about 95%
* [d5e3edb] four rounds of fuzzing of inputs and options: all crashes, hangs and memory errors
  fixed [bcdfbb0, 7494228]; embedded meshes and input values validated
* [6cc7990, 1c9093a] `test/fuzz/mmcfuzz.py`, an mmc fuzzer including the iMMC inputs, run on all four
  backends; the bugs it found were fixed [07c6fa5, 85f7ed9, 797cccf]

'''Packaging'''
* [73eb9f1] Linux releases built in manylinux2014 (glibc 2.17+); Octave 10 mex files for macOS/Windows
* [4e581d7] bundled zmat updated to 1.2.0; [b65562c] cJSON updated to 1.7.19

== # Introduction ==

Mesh-based Monte Carlo (MMC) traces photon packets through a tetrahedral mesh: at each
step a packet is attenuated by absorption, redirected by scattering, and reflected or
refracted when it crosses a face between elements of different refractive indices.
Accumulating the packets' weights gives a statistically exact solution of the
radiative transfer equation (RTE) [Fang2010]. Compared with voxel-based MC such as MCX:

* '''curved and oblique boundaries are exact''': reflection and refraction use the true
    surface normals, not voxel staircases, which matters at the skin, skull and other
    interfaces with a refractive-index change;
* '''thin and small structures''' (layers, vessels) are represented without making the
    whole grid finer, and the mesh can be coarse where the anatomy is simple;
* '''meshes can be shared with finite-element (FEM) solvers''' such as Redbird, so MC and
    diffusion results are directly comparable node by node.

The cost is that a mesh must be generated first (see Preparing a mesh)
and that tracing through tetrahedra is more work per step than through voxels.

MMC implements several variants of the method:

<pre>
mesh-based MMC
    What it does: traces and stores the fluence on the mesh (at nodes or elements)
        [Fang2010, Fang2012]
    Select with: `-M S` (`cfg.method='elem'`), or `P`, `H`, `B` on the CPU
dual-grid MMC (DMMC)
    What it does: traces on the mesh, stores the fluence on a fine voxel grid; 2-3x faster,
        and the default [Yan2019]
    Select with: `-M G` (`cfg.method='grid'`)
implicit MMC (iMMC)
    What it does: adds thin vessels, membranes and spheres given as radii/thicknesses on the
        mesh nodes, edges and faces, without meshing them [Yuan2021]
    Select with: `Mesh.MeshROI` (`cfg.noderoi`, `edgeroi`, `faceroi`) or
        `Mesh.MeshROIPrim`; CPU, OpenCL and CUDA
RT-MMC
    What it does: uses the ray-tracing cores of NVIDIA RTX GPUs through OptiX; also accepts
        surface-only meshes
    Select with: `-c optix` (`cfg.compute='optix'`)
wide-field MMC
    What it does: wide-field sources and detectors, by tessellating their domains into the
        mesh [Yao2016]
    Select with: done by MMCLAB/PMMC for wide-field sources
</pre>

'''The MCX family'''

<pre>
MCX (https://github.com/fangq/mcx)
    Description: voxel-based MC in CUDA (NVIDIA) or HIP (AMD); MCXLAB, PMCX and MCX Studio
MCX-CL (https://github.com/fangq/mcxcl)
    Description: the OpenCL version of MCX
MMC (https://github.com/fangq/mmc)
    Description: this package: mesh-based MC; MMCLAB and PMMC
iso2mesh (https://github.com/fangq/iso2mesh) / pyiso2mesh (https://pypi.org/project/iso2mesh/)
    Description: the mesh generators used to build MMC meshes from images and shapes
MCX Cloud (https://mcx.space/cloud)
    Description: run MCX and MMC in a web browser on cloud GPUs
Redbird (https://github.com/fangq/redbird) / Redbirdpy (https://github.com/fangq/redbirdpy)
    Description: diffuse optical tomography solvers that can use MMC as the forward model
</pre>

== # Choosing a backend ==

The same simulation runs on any of four backends, selected with `-c` on the command
line or `cfg.compute` in MMCLAB/PMMC:

<pre>
hardware
    `sse` (CPU): any x86_64 CPU (SSE4), and other CPUs
    `opencl`: any GPU or CPU with an OpenCL driver
    `cuda`: NVIDIA GPUs
    `optix`: NVIDIA GPUs; RTX GPUs use their ray-tracing cores
parallelism
    `sse` (CPU): OpenMP threads + SIMD
    `opencl`: GPU threads
    `cuda`: GPU threads
    `optix`: GPU threads
ray tracers (`-M`)
    `sse` (CPU): all: Plucker, Havel, Badouel, branchless Badouel, dual-grid
    `opencl`: branchless Badouel (`S`) and dual-grid (`G`)
    `cuda`: `S` and `G`
    `optix`: dual-grid (`G`) only
in the release packages
    `sse` (CPU): yes
    `opencl`: yes (default)
    `cuda`: mmc-optix packages
    `optix`: mmc-optix packages
in MMCLAB/PMMC releases
    `sse` (CPU): yes
    `opencl`: yes (default)
    `cuda`: build with `make cudamex/cudaoct`, or PMMC with `CMAKE_ARGS="-DBUILD_CUDA=ON"`
    `optix`: build with `make optixmex/optixoct`
</pre>

Feature support differs where noted below; everything else - all output types,
detectors, replay, adjoint Jacobians, RF, multiple sources, photon sharing - works on
every backend.

<pre>
implicit MMC (node/edge/face ROIs)
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: -
MCML-style weighting (`-m 1`)
    `sse`: yes
    `opencl`: -
    `cuda`: -
    `optix`: -
mesh (node/element) output
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: - (grid only)
`pencilarray`, `hyperboloid`, `ring` sources
    `sse`: yes
    `opencl`: -
    `cuda`: -
    `optix`: yes
photon replay (`-E file`, `-O J/L/P`)
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: -
photon sharing (`SrcNum` > 1)
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: -
detected photons
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: yes, without seeds, momentum or wide-field detectors
diffuse reflectance (`-X 1`)
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: -
per-node mua/mus' (`nodemua`, `nodemusp`)
    `sse`: yes
    `opencl`: yes
    `cuda`: yes
    `optix`: -
surface-only meshes (`MeshSurf`)
    `sse`: -
    `opencl`: -
    `cuda`: -
    `optix`: yes
identical photon paths for a seed
    `sse`: with one thread (`OMP_NUM_THREADS=1`)
    `opencl`: with `-A 0 -T 1`
    `cuda`: with `-A 0 -T 1`
    `optix`: -
</pre>

A backend that cannot run the requested settings stops with a message naming the
backend to use instead. Rules of thumb: use '''`opencl`''' by default, '''`sse`''' when no GPU
is available, '''`cuda`''' on NVIDIA GPUs when you build
MMC yourself, and '''`optix`''' on RTX GPUs for dual-grid forward simulations of large
meshes.

== # Installation and configuration ==

=== Getting MMC ===

<pre>
MMC (`mmc-<os>-<arch>-<version>.zip`)
    Contents: the `mmc` binary with the `sse` and `opencl` backends
    How to get it: https://mcx.space/wiki/?Get, https://github.com/fangq/mmc/releases
MMC-OptiX (`mmc-optix-<linux
    Contents: windows>-x64-<version>.zip`)
    How to get it: `mmc` with all four backends
    the same pages
MMCLAB (`mmclab-<os>-<arch>-<version>.zip`)
    Contents: MATLAB/Octave mex file (`sse`, `opencl`) and scripts
    How to get it: the same pages
MCX Studio all-in-one package
    Contents: mcx, mcxcl, mmc, the GUI and all MATLAB/Octave toolboxes
    How to get it: the same pages; on Windows also the installer or `winget install
        coti.mcxstudio`
PMMC (Python)
    Contents: `pip install pmmc` (`sse`, `opencl`)
    How to get it: https://pypi.org/project/pmmc/
nightly builds
    How to get it: https://mcx.space/nightly/
source code
    How to get it: `git clone https://github.com/fangq/mmc.git`
</pre>

=== The MCX Studio package ===

The all-in-one package unzips into one folder:

<pre>
MCXStudio/
├── mcxstudio            the GUI (MCXStudio.app on macOS)
├── mcxsuite_addpath.m   adds all the MATLAB/Octave toolboxes to the path
├── MATLAB/              mcxlab, mcxlabcl, mmclab, iso2mesh and redbird
└── MCXSuite/            mcx, mcxcl and mmc, each with its bin/ folder and helper scripts
</pre>

Add `MCXSuite/mmc/bin` (and `mcx/bin`, `mcxcl/bin`) to your `PATH`. In MATLAB or
Octave, `cd` into the `MCXStudio` folder, run `mcxsuite_addpath` (which also adds
the bundled iso2mesh), and `savepath` to keep
the paths for later sessions. On macOS, unzip the package under your home folder.

=== Requirements of each backend ===

* '''`sse`''': nothing beyond the operating system; the binary uses all CPU cores.
* '''`opencl`''': an OpenCL runtime, included in the NVIDIA, AMD and Intel graphics
    drivers and in macOS; on Linux, POCL (`pocl-opencl-icd`) adds the CPU as an OpenCL
    device. `mmc -L` lists the devices. See the MCX-CL README for details per vendor.
* '''`cuda`''' and '''`optix`''' (mmc-optix packages): an NVIDIA GPU with driver R515 or
    newer (CUDA 11.7). OptiX runs on Maxwell and newer GPUs, and uses hardware ray
    tracing on RTX GPUs. The CUDA toolkit and the OptiX SDK are only needed to build MMC.

The released Linux binaries run on glibc 2.17 or newer; the Octave mex files for
Windows and macOS need Octave 10 or newer.

=== Mesh tools ===

Building meshes needs [https://github.com/fangq/iso2mesh iso2mesh] in MATLAB/Octave
(included in the MCX Studio package, or `addpath` its folder; also `octave-iso2mesh` on Debian/Ubuntu/Fedora), or
[https://pypi.org/project/iso2mesh/ pyiso2mesh] in Python (`pip install iso2mesh`;
the PMMC examples need version 0.6.7 or newer). MMCLAB and PMMC call iso2mesh to
prepare the mesh for you.

=== Setting up MMCLAB and PMMC ===

<pre>
addpath('/path/to/mmclab');        % mmc.mex*, mmclab.m and the MMC helper functions
addpath('/path/to/iso2mesh');
info = mmclab('gpuinfo')           % list the OpenCL devices
</pre>

<pre>
python3 -m pip install pmmc pmcx iso2mesh jdata bjdata matplotlib
python3 -c "import pmmc; print(pmmc.gpuinfo())"
</pre>

`pmcx` provides the utility functions PMMC re-exports (`detweight`, `meanpath`, ...).

=== Building from source ===

<pre>
git clone https://github.com/fangq/mmc.git
cd mmc/src
make                        # mmc with sse + opencl, in ../bin
make cuda                   # mmc with sse + opencl + cuda (needs the CUDA toolkit)
make optix OPTIXHOME=/path/to/optix-sdk   # all four backends (needs CUDA 11+ and the OptiX 7.5 SDK)
make mex / make oct         # MMCLAB for MATLAB / Octave (sse + opencl), in ../mmclab
make cudamex / cudaoct      # MMCLAB with cuda
make optixmex / optixoct    # MMCLAB with cuda and optix
make unittest               # host unit tests
</pre>

`make sse`, `make omp` and `make release` build CPU-only variants, and
`make ... -f makefile_sfmt` uses the SFMT random number generator. CMake works too
(`cmake -S src -B build [-DBUILD_CUDA=ON]`), and is what PMMC uses:
`pip install .` in `pmmc/`, or `CMAKE_ARGS="-DBUILD_CUDA=ON" pip install .` for the
`cuda` backend. On macOS with clang, install OpenMP first (`brew install libomp`).
On Windows, MMC builds with MSYS2/MinGW, and the CUDA/OptiX backends with MSVC
(`.github/windows_release_builder/build_mmc_optix.sh`).

== # Quick start ==

'''Command line.''' Run a built-in benchmark, then a simulation from an input file:

<pre>
mmc -L                                   # list the OpenCL devices
mmc -Q dmmc-cube60 -n 1e6                # built-in benchmark (also dmmc-cube60b, edgeimmc, nodeimmc, faceimmc)
cd examples/onecube
mmc -f onecube.json -n 1e5               # mesh read from node_onecube.dat, elem_onecube.dat, ...
mmc -f onecube.json -n 1e5 -c sse        # the same on the CPU
mmc -f onecube.json --dumpjson onecube_all.json   # one self-contained JSON, mesh included
</pre>

'''MATLAB/Octave''' (with iso2mesh):

<pre>
cfg.nphoton = 1e6;
[cfg.node, face, cfg.elem] = meshabox([0 0 0], [60 60 30], 6);   % a 60x60x30 mm box mesh
cfg.elemprop = ones(size(cfg.elem, 1), 1);    % every element is medium 1
cfg.srcpos = [30 30 0];
cfg.srcdir = [0 0 1];
cfg.prop = [0 0 1 1; 0.005 1 0 1.37];          % [mua(1/mm), mus(1/mm), g, n] per label
cfg.tstart = 0; cfg.tend = 5e-9; cfg.tstep = 5e-10;
cfg.method = 'elem';                           % fluence at the mesh nodes
flux = mmclab(cfg);                            % flux.data: #nodes x #time gates
plotmesh([cfg.node log10(sum(flux.data, 2))], cfg.elem, 'y=30', 'edgecolor', 'none');
</pre>

'''Python''' (with pyiso2mesh):

<pre>
import numpy as np, iso2mesh as i2m, pmmc
cfg = {'nphoton': 1e6, 'srcpos': [30, 30, 0], 'srcdir': [0, 0, 1],
       'prop': [[0, 0, 1, 1], [0.005, 1, 0, 1.37]],
       'tstart': 0, 'tend': 5e-9, 'tstep': 5e-10, 'method': 'elem'}
cfg['node'], face, cfg['elem'] = i2m.meshabox([0, 0, 0], [60, 60, 30], 6)
cfg['elemprop'] = np.ones(cfg['elem'].shape[0])
res = pmmc.mmclab(cfg)                         # res['flux']: #nodes x #time gates
</pre>

'''GUI:''' start `mcxstudio`, open an MMC input file (a JSON file with a mesh), check the
mesh, source and detectors in the 3-D view, and press '''Run'''.

== # Preparing a mesh ==

An MMC mesh is a set of nodes (`node`, Nn x 3, in mm by default) and tetrahedra
(`elem`, Ne x 4, 1-based node indices), with a medium label for each element
(`elemprop`, Ne x 1, or a 5th column of `elem`). Labels start from 1; label 0 is the
space outside the mesh.

'''Making meshes with iso2mesh''' (MATLAB/Octave, or the same names in pyiso2mesh):

<pre>
simple shapes
    Function: `meshabox`, `meshasphere`, `meshacylinder`, `meshunitsphere`, and
        `mergemesh`/`surfboolean` to combine them
a segmented (label) volume, such as an MRI/CT segmentation
    Function: `v2m` / `vol2mesh` (labels become `elem(:,5)`)
closed surface meshes, such as from CAD or FreeSurfer
    Function: `s2m` / `surf2mesh`
</pre>

The element size sets the spatial resolution of mesh-based output; with dual-grid MMC
(the default), the output resolution is set by the grid instead, so the mesh only
needs to follow the boundaries. [https://github.com/COTILab/BlenderPhotonics BlenderPhotonics]
builds meshes interactively in Blender and runs MMC on them.

'''Requirements MMC checks for'''

* '''consistent orientation''': all tetrahedra must have the same node ordering.
    MMCLAB/PMMC reorient them (`meshreorient`); for the command line, save meshes with
    `savemmcmesh.m`, which does the same.
* '''no degenerate elements''': elements of zero volume are rejected.
* '''the source inside the mesh''': MMC needs the element that contains the source
    (`cfg.e0`, `Mesh.InitElem`); MMCLAB/PMMC find it with `tsearchn`. For a source
    outside the mesh, set `cfg.e0` to `'>'` to move the source along `srcdir` to the
    first surface it meets (`'<'` backwards, `'-'` both).

'''Wide-field sources''' (`planar`, `pattern`, `fourier`, ...) must start ''outside'' the
mesh. MMCLAB/PMMC then add the source's footprint to the mesh as extra elements with
negative labels (`mmcsrcdomain`, `mmcaddsrc`), from which photons are launched;
`mmcadddet` does the same for wide-field detectors. If the source sits on the mesh
surface, mmclab stops with an explanation, because the GPU would otherwise launch every
photon from one element.

'''Using a mesh on the command line.''' Either save it as files and give their name stub
in `Mesh.MeshID`:

<pre>
savemmcmesh('mymesh', node, elem);   % node_mymesh.dat, elem_mymesh.dat, facenb_..., velem_...
</pre>

or embed it in the JSON input (`Mesh.MeshNode`, `Mesh.MeshElem`, see
JSON input file format). The easiest way to get a complete
input file is to set the simulation up in MATLAB and write it out with
`mmc2json(cfg, 'input.json')`, or to run `mmc ... --dumpjson`.

'''Surface-only meshes (OptiX).''' `cfg.face` (`Mesh.MeshSurf`): triangles given as
`[n1 n2 n3 inner_label outer_label]`, with normals pointing from the inner to the outer
medium, replace the tetrahedral mesh when `cfg.compute='optix'`; no volume mesh is
needed. Example: `demo_optix_mmc.m`.

== # Key concepts ==

'''Units and positions.''' Node coordinates, source and detector positions and radii are
in the length unit of the mesh, mm by default; set `Mesh.LengthUnit`/`-u`/
`cfg.unitinmm` if the mesh is in another unit (e.g. 10 for cm). Absorption (mua) and
scattering (mus) coefficients are always in 1/mm.

'''Media.''' The optical properties of label i are the (i+1)-th row of `cfg.prop`
(`Domain.Media`): `[mua, mus, g, n]`, with mus the scattering (not reduced scattering)
coefficient. The first row is label 0, the space outside the mesh; its refractive index
can also be set with `cfg.nout`. Per-node mua/mus' (`cfg.nodemua`, `cfg.nodemusp`) can
replace the per-label values, as image reconstruction needs.

'''Where the fluence is stored.''' With the default dual-grid method (`-M G`,
`cfg.method='grid'`), the output is a voxel volume covering the bounding box of the
mesh, with voxels of size `--gridsize`/`cfg.steps` (1 mm by default). With the
mesh tracers (`-M S`, `cfg.method='elem'`, ...), the output is stored on the mesh: at the
nodes (`-C 1`, `cfg.basisorder=1`, the default) or in the elements (`-C 0`). Either can
be plotted with iso2mesh (`plotmesh`, `qmeshcut`) or any volume viewer.

'''Time gates, photons and normalization''' work as in MCX: the time window [T0, T1] is
divided into gates of width Dt, results are averages over the simulated photons (noise
falls as 1/sqrt(photons)), and the output is normalized to a unit source by default
(the Green's function; `-U 0` for raw weights).

'''Boundaries.''' `-b`/`cfg.isreflect` sets the condition at the mesh surface: 1 Fresnel
reflection (the default, using `nout` outside), 0 index-matched (photons leave), 2 total
absorption, 3 mirror reflection. Interior interfaces between media of different
refractive indices always reflect and refract.

== # Running the mmc command line ==

<pre>
mmc -f input.json [options]
mmc -Q dmmc-cube60 -n 1e7
</pre>

`mmc` (or `mmc -h`) prints all options; the most useful ones are below.

<pre>
`-f file`, `--input`
    Meaning: read a JSON input file (inline JSON if the string starts with `{`, stdin with
        `-f -`)
`-Q name`, `--bench`
    Meaning: run a built-in benchmark: dmmc-cube60, dmmc-cube60b, edgeimmc, nodeimmc,
        faceimmc
`-N name`, `--net`
    Meaning: download and run a simulation shared on NeuroJSON.io
        (https://neurojson.io/db/mmc)
`-j '{...}'`, `--json`
    Meaning: modify any setting of the input with a JSON fragment
`-n num`, `--photon`
    Meaning: number of photons
`-E seed`, `--seed`
    Meaning: RNG seed, or a `.mch` file to replay its detected photons
`-P n`, `--replaydet`
    Meaning: replay only the photons of detector n
`-c name`, `--compute`
    Meaning: backend: `sse`, `opencl` (default), `cuda`, `optix`
`-G n` or `-G '0110'`, `--gpu`
    Meaning: device n (from 1, see `-L`), several devices as a 0/1 mask, or `-1` for the CPU
        (`sse`)
`-W '50,50'`, `--workload`
    Meaning: relative share of photons for each device
`-L`, `--listgpu`
    Meaning: list the OpenCL devices
`-M G`, `--method`
    Meaning: ray tracer: `G` dual-grid (default), `S` branchless Badouel, and on the CPU `P`
        Plucker, `H` Havel, `B` Badouel
`--gridsize mm`
    Meaning: voxel size of the dual-grid output
`-C 1/0`, `--basisorder`
    Meaning: mesh output at the nodes (1) or in the elements (0)
`-u mm`, `--unitinmm`
    Meaning: length unit of the mesh in mm
`-b 1/0`, `--reflect`
    Meaning: Fresnel reflection at the mesh surface (1) or none (0)
`-U 1/0`, `--normalize`
    Meaning: normalize the output (1) or save raw weights (0)
`-m 0/1`, `--mc`
    Meaning: MCX-style (0) or MCML-style (1, CPU only) photon weighting
`-e w`, `--minenergy`
    Meaning: Russian roulette threshold (default 1e-6)
`-V 0/1`, `--specular`
    Meaning: include the specular reflection for a source outside the mesh
`-k 1/0`, `--voidtime`
    Meaning: for a source outside the mesh, start the clock at launch (1) or at entry (0)
`--srcid n`
    Meaning: with several sources: 0 all together, -1 each separately, -2 also run detectors
        as sources, n only source n
`-s name`, `--session`
    Meaning: name of the output files
`-O X`, `--outputtype`
    Meaning: output quantity, see Output quantities
`-F fmt`, `--outputformat`
    Meaning: `ascii` (default, `.dat`), `bin`, `mc2`, `nii`, `hdr`, `jnii`, `bnii`
`-Z zlib`, `--zip`
    Meaning: compression for jnii/bnii
`-d 1`, `--savedet`
    Meaning: save detected photons (off by default in mmc)
`-x 1`, `--saveexit`
    Meaning: also save exit positions and directions; 2: a wide-field detector image
`--momentum 1`
    Meaning: also save momentum transfer (for DCS)
`-X 1`, `--saveref`
    Meaning: save the diffuse reflectance on the mesh surface
`-q 1`, `--saveseed`
    Meaning: save the RNG seeds of detected photons, for replay
`-H num`, `--maxdetphoton`
    Meaning: the most detected photons kept
`-S 1/0`, `--save2pt`
    Meaning: save the volume output (1) or skip it (0)
`--dumpjson [file]`
    Meaning: print (or save) the complete input, mesh included, as one JSON file
`-A 1/0`, `--autopilot`
    Meaning: let mmc choose the GPU thread count, with fast math on OpenCL (1), or use
        `-t`/`-T` (0)
`--atomic 1/0`, `--kernel file.cl`
    Meaning: atomic accumulation; a custom OpenCL kernel
`-D flags`, `--debug`
    Meaning: debug output: `P` progress bar, `T` timing, `M` trajectories, and ray-tracing
        details (`S`, `C`, `B`, `W`, `D`, `I`, `O`, `X`, `A`, `R`, `E`)
`--debugphoton n`
    Meaning: print the debug output for photon n only
</pre>

'''Recipes'''

<pre>
mmc -f input.json -c opencl -G 1                     # first OpenCL device
mmc -f input.json -c sse                             # all CPU cores
mmc -f input.json -c cuda                            # NVIDIA GPU (mmc-optix build)
mmc -f input.json -c optix -M G --gridsize 0.5       # RT-MMC, 0.5 mm output voxels
mmc -f input.json -M S -C 1 -F jnii                  # fluence at the mesh nodes, as JNIfTI
mmc -f input.json -d 1 -q 1 -s base                  # save detected photons and seeds (base.mch)
mmc -f input.json -E base.mch -O J -P 1 -s jac       # replay them: mua Jacobian of detector 1
</pre>

=== Debug output ===

With the ray-tracing debug flags (`-D S` photon moves, `A` accumulation, `E` exits, and
`C`, `B`, `W`, `D`, `I`, `O`, `X`, `R`), mmc prints one JSON object per line:

<pre>
{"id":0, "pre-add": {"p":[2.000000, 8.000000, 0.407062], "lmin":1.333540e+00, "eid":3, "l":8.814591e-01}}
{"id":0, "terminate": {"p":[2.152560, 10.000000, 1.270147], "v":[-0.623828, 0.771854, -0.122802], "w":0.981131, "eid":-6}}
</pre>

`id` is the photon index and the key names the event (`move`, `hit`, `enter`, `leave`,
`exit-boundary`, `add`, `roulette-die`, `terminate-by-time`, the implicit-MMC
`implicit-node`/`edge`/`face`, ...); the record holds the position `p`, direction `v`,
weight `w`, element `eid` and step length `l` where they apply. Use it on the CPU
(`-c sse`) with a few photons, or with `--debugphoton n` to follow one photon. To plot
photon paths, save trajectories instead (`-D M`, or the `traj` output of MMCLAB).

== # JSON input file format ==

An MMC input file has the same structure as an MCX input, with a `Mesh` section in
place of the voxel grid (`Shapes` is accepted as its name too):

<pre>
{
    "Mesh": {
        "MeshID": "onecube",
        "InitElem": 3
    },
    "Session": {"Photons": 100000, "RNGSeed": 17182818, "ID": "onecube"},
    "Forward": {"T0": 0.0, "T1": 5.0e-09, "Dt": 5.0e-10},
    "Optode": {
        "Source": {"Pos": [2.0, 8.0, 0.0], "Dir": [0.0, 0.0, 1.0]},
        "Detector": [{"Pos": [2.0, 6.0, 0.0], "R": 1.0}]
    }
}
</pre>

'''`Mesh`'''

<pre>
`MeshID`
    Meaning: name stub of the mesh files `node_<id>.dat`, `elem_<id>.dat`,
        `facenb_<id>.dat`, `velem_<id>.dat` (and `prop_<id>.dat` when `Domain.Media` is
        absent), written by `savemmcmesh.m`
`MeshNode`
    Meaning: the nodes, Nn x 3, as a JSON array or a JData-encoded `single` array
`MeshElem`
    Meaning: the elements, Ne x 5: 4 node indices and the label (`int32`)
`MeshSurf`
    Meaning: OptiX surface-only mesh, Nf x 5: 3 node indices, inner and outer labels
`MeshROI`
    Meaning: implicit MMC: node radii (Nn x 1), edge radii (Ne x 6) or face thicknesses (Ne
        x 4)
`MeshROIPrim`
    Meaning: implicit MMC: capsules and spheres given by coordinates, N x 7
        [Ax,Ay,Az,Bx,By,Bz,r] (a sphere if B equals A), not tied to the mesh
`InitElem`
    Meaning: the element containing the source
`LengthUnit`
    Meaning: length unit of the mesh in mm
</pre>

'''`Domain`''': `Media` (array of `{"mua", "mus", "g", "n"}`, the first for label 0),
`LengthUnit`, and `Step` (the dual-grid voxel size).

'''`Session`''': `ID`, `Photons`, `RNGSeed`, `OutputType`, `OutputFormat`, `RayTracer`
(`-M` letter), `BasisOrder`, `DoMismatch`, `DoNormalize`, `DoSaveVolume`,
`DoPartialPath`, `DoSaveExit`, `DoSaveSeed`, `DoDCS`, `DoSpecular`, `DebugFlag`.

'''`Forward`''': `T0`, `T1`, `Dt` (seconds), `N0` (refractive index outside the mesh) and
`Omega` (modulation angular frequency for RF outputs).

'''`Optode`''': `Source` (`Type`, `Pos`, `Dir`, `Param1`, `Param2`, `Pattern`, `SrcNum`;
or an array of sources) and `Detector` (an array of `{"Pos", "R"}`, plus `"Dir"` for the
adjoint outputs; `Param1`/`Param2` for a wide-field detector).

`mmc2json(cfg, 'file.json')` in MATLAB converts an MMCLAB `cfg`, mesh included, into an
input file, and `mmc --dumpjson` writes the complete input of any run.

== # Sources ==

MMC supports the same source types as MCX (see the MCX README for their parameters):
`pencil`, `isotropic`, `cone`, `gaussian`, `planar`, `pattern`, `fourier`, `arcsine`,
`disk`, `fourierx`, `fourierx2d`, `zgaussian`, `line`, `slit`, `pencilarray`, `pattern3d`,
`hyperboloid` and `ring`. Positions and sizes are in mesh units. The OpenCL and CUDA
backends do not implement `pencilarray`, `hyperboloid` and `ring`; use `sse` or `optix`
for them.

* '''Point and pencil sources''' must lie inside the mesh (in element `InitElem`/`cfg.e0`),
    or outside with `cfg.e0='>'` so MMCLAB moves them to the surface.
* '''Wide-field sources''' (`planar`, `pattern`, `fourier`, `fourierx`, `fourierx2d`,
    `pattern3d`, ...) must start outside the mesh; MMCLAB/PMMC attach their footprint
    to the mesh (see Preparing a mesh). Example:
    `demo_sfdi_2layer.m`.
* '''Several sources''': one row per source in `cfg.srcpos`/`cfg.srcdir`, or an array of
    `Optode.Source` objects; `cfg.srcid` chooses whether they are simulated together or
    separately. Example: `demo_mmclab_multisrc.m`.
* '''Photon sharing''': a `pattern` source with `cfg.srcnum` > 1 simulates many patterns
    in one run. Example: `demo_photon_sharing.m`.

== # Detectors and detected photon data ==

Detectors are disks on the mesh surface, `cfg.detpos = [x y z radius]`
(`Optode.Detector`). On the command line, detected photons are saved only with `-d 1`.
For each detected photon, MMC saves the detector ID, the scattering count and partial
path length in each medium, the momentum transfer (`--momentum 1`, `cfg.ismomentum=1`),
the exit position and direction (`-x 1`, `cfg.issaveexit=1`) and the initial weight.
The helper functions `mmcdetweight`, `mmcmeanpath`, `mmcmeanscat`, `mmcdettpsf` and
`mmcdettime` (MATLAB), or `pmmc.detweight`, `pmmc.meanpath`, ... (Python), turn them
into detector readings for any optical properties. In the raw array (`detp.data` in
MMCLAB, or the `.mch` file read by `loadmch`), each photon is a row of: detector ID,
the scattering counts and partial path lengths in each medium, the momentum transfers
(with `--momentum 1`), the exit position and direction (with `-x 1`), and the initial
weight.

A '''wide-field (camera) detector''' is a quadrilateral given by `cfg.detpos` and two
edge vectors with pixel counts, `cfg.detparam1 = [ux uy uz Nx]` and
`cfg.detparam2 = [vx vy vz Ny]`; with `-x 2` (`cfg.issaveexit=2`) MMC returns the
image of the exiting weight instead of per-photon records. Example: `demo_wide_det.m`.

== # Output quantities, files and formats ==

`-O` (`Session.OutputType`, `cfg.outputtype`) selects what the volume output holds:

<pre>
`X`
    Name: `flux`
    Quantity: fluence rate, 1/(mm^2 s) per unit source (default)
`F`
    Name: `fluence`
    Quantity: fluence, integrated over each time gate
`E`
    Name: `energy`
    Quantity: energy deposited
`J`
    Name: `jacobian`
    Quantity: replay: mua Jacobian
`L`
    Name: `wl`
    Quantity: replay: weighted path length (mua Jacobian)
`P`
    Name: `wp`
    Quantity: replay: weighted scattering count (mus Jacobian)
`R`
    Name: `rf`
    Quantity: with a modulation frequency, no replay: complex fluence rate
`A`, `D`, `U`, `V`
    Name: `adjoint`, `adjoint_dcoeff`, `adjoint_mus`, `adjoint_musp`
    Quantity: adjoint Jacobians
`W`, `Q`
    Name: `adjoint_mua_d`, `adjoint_mua_musp`
    Quantity: two adjoint Jacobians together
</pre>

'''Files of the command line''': the volume output `<session>.<ext>` in the format set by
`-F` (`.dat` text by default, or `.bin`, `.mc2`, `.nii`, `.jnii`, `.bnii`; RF runs also
write the imaginary part to `<session>_imag.<ext>`), the detected photons
`<session>.mch` (or `<session>_detp.jdt`/`.jdb` with `-F jnii`/`bnii`), the surface
reflectance `<session>_dref.<ext>` with `-X 1`, and trajectories with `-D M`. With
mesh output, each time gate is a block of values at the nodes (or elements); with
dual-grid output, it is a voxel volume.

The diffuse reflectance (`-X 1`) is stored for each triangle of the mesh surface, in
the order of `faceneighbors(elem, 'rowmajor')` (iso2mesh), as `#triangles x #gates`.

=== Units of the output ===

With the default normalization (`-U 1`), the output is the response to a source that
emits one unit of light (one joule, or one photon) in an impulse at t=0, so each
quantity can be read in energy or in photon units:

<pre>
fluence rate (`flux`, default)
    Energy units: W/mm^2 per J emitted
    Photon units: 1/(mm^2 s) per photon emitted
fluence (`fluence`)
    Energy units: J/mm^2 per J emitted
    Photon units: 1/mm^2 per photon emitted
energy (`energy`)
    Energy units: J absorbed in the node, element or voxel per J emitted
    Photon units: the fraction of the photons absorbed there
</pre>

To obtain absolute values, multiply by the energy (or number of photons) your source
delivers; for a continuous source of power P, the CW fluence rate is P times the
time-integrated (summed over gates and multiplied by Dt) fluence rate. Divide the
`energy` output by the volume it represents (the nodal or element volume, or the voxel volume) for the absorbed energy density;
`mcxfluence2energy.m` converts a fluence output into energy deposition. See [Fang2010] and Lihong Wang's *Biomedical Optics: Principles
and Imaging*, chapter 5, for the definitions.

'''Fluence rate and current density.''' The fluence rate is a scalar: the light passing
through a small sphere around a point from all directions. The current density J is a
vector: the net light flowing through a surface in a given direction; it has the same
units. What a detector or the diffuse reflectance measures is the current density
flowing out of the mesh: the total weight of the detected photons divided by the
detector area gives J at the detector, as a function of time when the photons are
binned by their time of flight, and `-X 1` records it on every surface triangle.

== # Reading MMC data in MATLAB and Python ==

<pre>
flux = load('onecube.dat');                 % -F ascii: one value per node (or element) per line
nii  = loadjd('onecube.jnii');              % JSONLab: .jnii, .bnii, .jdt, .jdb
[data, header] = loadmch('onecube.mch');    % detected photons (mmc/matlab)
[node, elem] = readmmcmesh('onecube');      % read node_onecube.dat and elem_onecube.dat
</pre>

Mesh output is plotted with iso2mesh: `plotmesh([node log10(fluence)], elem, 'y=30')`
shows a cross section, `plotmesh([node fluence], elem, 'x<30 & y>30')` a cut-away
view, and `qmeshcut` returns the values on any plane for your own plots. Dual-grid
output is a regular volume, shown with `imagesc`, `islicer` or any volume viewer.

<pre>
import jdata as jd, pmmc
nii  = jd.load('onecube.jnii')              # needs bjdata for .bnii/.jdb
data, header = pmmc.loadmch('onecube.mch')[:2]
</pre>

== # Using MMCLAB in MATLAB and Octave ==

MMCLAB runs MMC inside MATLAB or GNU Octave, and does the mesh bookkeeping for you:
it reorients the elements, computes face neighbors and element volumes, finds the
source element, and attaches wide-field source domains, before calling the `mmc` mex
file.

<pre>
[flux, detp, ncfg, seeds, traj] = mmclab(cfg);
ncfg = mmclab(cfg, 'prep');            % only preprocess, to reuse the result
flux = mmclab(cfg, 'sse');             % run on a given backend: 'sse', 'opencl', 'cuda', 'optix'
info = mmclab('gpuinfo');
</pre>

`help mmclab` documents every `cfg` field, its default and the outputs.

<pre>
`flux.data`
    Content: dual-grid: Nx x Ny x Nz x Ng voxels; mesh output: #nodes (or #elements) x Ng
`flux.dref`
    Content: surface diffuse reflectance, #triangles x Ng, with `cfg.issaveref=1`
`flux.jmua`, `.jd`, `.jmus`, `.jmusp`
    Content: adjoint Jacobians
`detp`
    Content: detected photons: `detid`, `nscat`, `ppath`, `mom`, `p`, `v`, `w0`
`ncfg`
    Content: the preprocessed cfg, to pass to later runs
`seeds`, `traj`
    Content: RNG seeds of detected photons; trajectories (`id`, `pos`, `eid`)
</pre>

'''Helper functions''' (`mmc/matlab`): mesh I/O (`savemmcmesh`, `readmmcmesh`,
`readmmcnode`, `readmmcelem`, `readmmcface`), sources and detectors (`mmcsrcdomain`,
`mmcaddsrc`, `mmcadddet`, `mmcraytrace`), detected photons (`mmcdetweight`,
`mmcmeanpath`, `mmcmeanscat`, `mmcdettpsf`, `mmcdettime`, `loadmch`), Jacobians
(`mmcjacobian`, `mmcjmua`, `mmcjmus`), regular meshes (`genT5mesh`, `genT6mesh`), and the
analytical diffusion solutions of spheres (`sphdiffusion*`); `mmc2json` is in
`mmclab`.

'''Examples''' in `mmclab/example`: `demo_mmclab_basic.m` (start here),
`demo_example_onecube.m`, `demo_head_atlas.m` (a brain atlas mesh),
`demo_dmmc_sphshells.m` (dual-grid MMC), `demo_immc_basic.m`, `demo_immc_vessel.m`
(implicit MMC), `demo_sfdi_2layer.m`, `demo_wide_det.m` (wide-field sources and
detectors), `demo_example_replay.m`, `demo_mmclab_mesh_adjoint.m`,
`demo_mmclab_mesh_rf_adjoint.m` (Jacobians), `demo_photon_sharing.m`,
`demo_mmclab_multisrc.m`, `demo_optix_mmc.m` (RT-MMC) and `demo_compare_mmc_mcx.m`.

== # Using PMMC in Python ==

PMMC runs MMC from Python:

<pre>
import pmmc
res = pmmc.run(cfg)          # cfg: a dict with the MMCLAB field names, mesh already prepared
res = pmmc.mmclab(cfg)       # with MMCLAB's preprocessing (needs pyiso2mesh)
print(pmmc.gpuinfo())
</pre>

`pmmc.mmclab` reorients the mesh, finds the source element and attaches wide-field
source domains as `mmclab.m` does; `pmmc.mmcsrcdomain`, `mmcaddsrc`, `mmcadddet` and
`mmcraytrace` are also available, together with the detected-photon and plotting
utilities from `pmcx`. All MMCLAB examples have Python versions in `pmmc/example`.
The wheels include the `sse` and `opencl` backends; `CMAKE_ARGS="-DBUILD_CUDA=ON"
pip install .` in `pmmc/` adds `cuda`.

== # Using MCX Studio ==

MCX Studio 2, included in the MCX Studio packages, is the graphical front end for MCX,
MCX-CL and MMC (source code in MCX's `mcxstudio2/` folder). It opens MMC input files,
shows the mesh surface, sources, detectors, the fluence and photon paths in its 3-D
view, and runs mmc when the domain is a mesh; the `cpu` processor runs mmc's `sse`
backend:

<pre>
mcxstudio onecube.json
xvfb-run -a mcxstudio --run --quit --processor cpu onecube.json   # headless
</pre>

== # Advanced features ==

=== Photon replay and Jacobians ===

Save the seeds of the detected photons, then replay them to build Jacobians
[Yao2018]. In MMCLAB:

<pre>
[flux, detp, ncfg, seeds] = mmclab(cfg);
ncfg.seed = seeds.data;
ncfg.outputtype = 'jacobian';      % or 'wl', 'wp'
ncfg.replaydet = 1;                % detector 1 (0: all detectors summed)
jac = mmclab(ncfg);
</pre>

On the command line: `mmc -f in.json -d 1 -q 1 -s base`, then
`mmc -f in.json -E base.mch -O J -P 1`. Example: `demo_example_replay.m`.

=== Adjoint Jacobians ===

Give each detector a direction (`cfg.detdir`, `Optode.Detector.Dir`) and use output
type `adjoint`, `adjoint_dcoeff`, `adjoint_mus`, `adjoint_musp`, `adjoint_mua_d` or
`adjoint_mua_musp`: every detector is also run as a source, and the Jacobians of all
source-detector pairs are formed from the forward and adjoint fields, in CW or the
frequency domain. On the mesh, `cfg.adjointmode` chooses the full finite-element form
(0, the default) or the nodal approximation (1). Examples:
`demo_mmclab_mesh_adjoint.m`, `demo_mmclab_mesh_rf_adjoint.m`. Redbird uses these
outputs, with `cfg.nodemua`/`cfg.nodemusp`, for image reconstruction.

=== Frequency-domain (RF) simulations ===

With a modulation frequency (`cfg.omega` in rad/s, `Forward.Omega`), MMC tracks complex
photon weights and returns the complex fluence rate, on every backend.

=== Dual-grid MMC ===

The default method traces photons on the mesh but stores the fluence on a voxel grid,
which is faster and gives a uniform output resolution; set the grid with `cfg.steps` or
`--gridsize`. `demo_dualmesh_output.m` shows how to map the result back to the mesh.

=== Implicit MMC ===

Thin vessels, membranes and small spheres can be added without meshing them: give a
radius per node (`cfg.noderoi`), per element edge (`cfg.edgeroi`, Ne x 6) or a thickness
per element face (`cfg.faceroi`, Ne x 4) [Yuan2021]. Implicit MMC runs on the CPU
(`sse`), OpenCL and CUDA backends, not on OptiX. Examples: `demo_immc_basic.m`, `demo_immc_vessel.m`; benchmarks `edgeimmc`,
`nodeimmc`, `faceimmc`.

=== Comparing with MCX ===

`demo_compare_mmc_mcx.m` runs the same problem with MMC and MCX; meshes and voxel
volumes of the same domain give statistically matching results, which the test suite
checks.

== # Performance and best practices ==

* '''Use a GPU backend for large runs''' (`opencl`, `cuda` or `optix`); keep `sse` for
    MCML-style weighting and runs without a GPU.
* '''Prefer dual-grid output''' (the default) when you need a fluence map: it is faster
    than mesh output and its resolution does not depend on the mesh.
* '''Mesh only as finely as the geometry needs.''' Very small or badly shaped elements
    slow down tracing without improving the result.
* '''Preprocess once.''' `mmclab(cfg, 'prep')` returns a cfg with the mesh data computed;
    reuse it for repeated runs.
* '''Reproducing a run''': on the CPU, use one thread (`OMP_NUM_THREADS=1`) and the same
    seed; on the GPU, `-A 0 -T 1` gives the same photon paths.
* '''More photons, less noise''': noise falls as 1/sqrt(photons).

== # Testing ==

<pre>
host unit tests
    What it checks: parsers, mesh code, ray tracers, and the CUDA/OptiX backends where built
    How to run: `make unittest` in `src`
functional tests
    What it checks: absorbed fraction, detected photons and outputs of the examples, on each
        backend that is built
    How to run: `test/testmmc.sh`
built-in examples
    What it checks: example cases exported from the MMCLAB demos, against golden values
    How to run: `test/testbuiltin.sh`
trajectory tests
    What it checks: step-by-step photon trajectories on the CPU and OpenCL, against golden
        data and MCX
    How to run: `test/testtraj.sh`
PMMC tests
    What it checks: the Python module
    How to run: `python3 -m unittest discover -s test` in `pmmc/`
</pre>

== # Known limitations ==

* MMC supports linear tetrahedral elements only.
* Optical properties are defined per element (per label), except the per-node
    mua/mus' (`cfg.nodemua`, `cfg.nodemusp`) used for image reconstruction.
* Implicit MMC (node/edge/face ROIs) does not run on OptiX, MCML-style weighting runs
    on the CPU (`sse`) only, and the OptiX backend has the restrictions listed in
    Choosing a backend.
* MMC implements the mesh-based Monte Carlo method as described in the papers below,
    not every possible variant of it. If you compare MMC with other software in a
    publication, please contact the author to make sure the details of the
    implementation are understood correctly.

== # Getting help ==

* documentation, FAQ and tutorials: https://mcx.space/wiki/
* `help mmclab`, `pmmc/README.md`, and the examples
* questions and discussion: the [https://groups.google.com/g/mmc-users mmc-users mailing list] for MMC,
    and the [https://groups.google.com/g/mcx-users mcx-users mailing list] for all MCX software
* bug reports and feature requests: https://github.com/fangq/mmc/issues
* release announcements: the [https://groups.google.com/g/mcx-news mcx-news mailing list]

== # Contributing ==

MMC is open-source software under the GPLv3 license: you may use, modify and share
it, as long as derived work is also released under the GPLv3. Bug fixes and new
features are welcome:

* fork https://github.com/fangq/mmc, make your change in your fork, keep it small
    and local, and open a pull request; describe what it fixes and how you tested it;
* or send the change as a `git diff` to the mailing list or a GitHub issue.

Please format C/C++ code with `make pretty` (astyle) in `src`, MATLAB code with
`mh_style --fix` (miss_hit) and Python code with `black` before submitting, and add a
test when you can. Contributions are acknowledged in the What's New section and the
ChangeLog.

== # Publications and how to cite ==

If you use MMC in your research, please cite [Fang2010], which introduced MMC,
[Fang2019], which describes the GPU-accelerated MMC, and the papers of the specific
features you use. The complete list is maintained at https://mcx.space/#publications.

* [Fang2010] Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing in
  Plücker coordinates," Biomed. Opt. Express, 1(1), 165-175 (2010).
  https://opg.optica.org/boe/abstract.cfm?uri=boe-1-1-165
* [Fang2019] Qianqian Fang and Shijie Yan, "Graphics processing unit-accelerated mesh-based
  Monte Carlo photon transport simulations," J. Biomed. Opt. 24(11), 115002 (2019).
  https://doi.org/10.1117/1.JBO.24.11.115002
* [Fang2012] Qianqian Fang and David R. Kaeli, "Accelerating mesh-based Monte Carlo method on
  modern CPU architectures," Biomed. Opt. Express, 3(12), 3223-3230 (2012).
  https://opg.optica.org/boe/abstract.cfm?uri=boe-3-12-3223
* [Yan2019] Shijie Yan, Anh Phong Tran, Qianqian Fang, "A dual-grid mesh-based Monte Carlo
  algorithm for efficient photon transport simulations in complex 3D media," J. Biomed. Opt.
  24(2), 020503 (2019). https://doi.org/10.1117/1.JBO.24.2.020503
* [Yuan2021] Yaoshen Yuan, Shijie Yan, Qianqian Fang, "Light transport modeling in highly
  complex tissues using the implicit mesh-based Monte Carlo algorithm," Biomed. Opt.
  Express, 12(1), 147-161 (2021). https://opg.optica.org/boe/fulltext.cfm?uri=boe-12-1-147
* [Yao2016] Ruoyang Yao, Xavier Intes, Qianqian Fang, "Generalized mesh-based Monte Carlo for
  wide-field illumination and detection via mesh retessellation," Biomed. Opt. Express,
  7(1), 171-184 (2016). https://opg.optica.org/boe/abstract.cfm?uri=boe-7-1-171
* [Yao2018] Ruoyang Yao, Xavier Intes, Qianqian Fang, "A direct approach to compute Jacobians
  for diffuse optical tomography using perturbation Monte Carlo-based photon 'replay',"
  Biomed. Opt. Express, 9(10), 4588-4603 (2018).
  https://opg.optica.org/boe/abstract.cfm?uri=boe-9-10-4588
* Qianqian Fang, "Comment on 'A study on tetrahedron-based inhomogeneous Monte-Carlo optical
  simulation'," Biomed. Opt. Express, 2(5), 1258-1264 (2011).
  https://opg.optica.org/boe/abstract.cfm?uri=boe-2-5-1258
* Yuxuan Zhang and Qianqian Fang, "BlenderPhotonics - an integrated open-source software
  environment for three-dimensional meshing and photon simulations in complex tissues,"
  J. Biomed. Opt. 27(8), 083014 (2022). https://doi.org/10.1117/1.JBO.27.8.083014
* Anh Phong Tran, Shijie Yan, Qianqian Fang, "Improving model-based fNIRS analysis using
  mesh-based anatomical and light-transport models," Neurophotonics, 7(1), 015008 (2020).
  https://doi.org/10.1117/1.NPh.7.1.015008
* Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon Migration in 3D Turbid
  Media Accelerated by Graphics Processing Units," Opt. Express, 17(22), 20178-20190 (2009),
  the MCX paper.

== # Acknowledgement ==

The MCX project, including MMC, is funded by the US National Institutes of Health (NIH) /
National Institute of General Medical Sciences (NIGMS) under grant R01-GM114365. We thank
the users and contributors who reported bugs and contributed code, credited in the What's
New section and the ChangeLog.

MMC uses the following open-source libraries:

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
* Comment:
 All the cruft you find here is public domain.  You don't have to
 credit anyone to use this code, but my personal request is that you mention
 Igor Pavlov for his hard, high quality work.

=== Miniz compression library ===

* Files: src/zmat/miniz/*
* Copyright 2013-2014 RAD Game Tools and Valve Software
* Copyright 2010-2014 Rich Geldreich and Tenacious Software LLC
* URL: https://github.com/richgel999/miniz
* License: MIT-license, https://github.com/richgel999/miniz/blob/master/LICENSE

=== SSE Math library by Julien Pommier ===

Copyright (C) 2007 Julien Pommier

This software is provided 'as-is', without any express or implied warranty. In
no event will the authors be held liable for any damages arising from the use
of this software.

Permission is granted to anyone to use this software for any purpose, including
commercial applications, and to alter it and redistribute it freely, subject to
the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim
that you wrote the original software. If you use this software in a product, an
acknowledgment in the product documentation would be appreciated but is not
required.

2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.

(this is the zlib license)

=== cJSON library by Dave Gamble ===

Copyright (c) 2009 Dave Gamble

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the “Software”), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=== SFMT library by Mutsuo Saito, Makoto Matsumoto and Hiroshima University ===

Copyright (c) 2006,2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
University. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

* Neither the name of the Hiroshima University nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=== drand48_r port for libgw32c by Free Software Foundation ===

Copyright (C) 1995, 1997, 2001 Free Software Foundation, Inc. This file is part
of the GNU C Library. Contributed by Ulrich Drepper
&lt;drepper@gnu.ai.mit.edu&gt;, August 1995.

The GNU C Library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU C Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
more details.

You should have received a copy of the GNU Lesser General Public License along
with the GNU C Library; if not, write to the Free Software Foundation, Inc., 59
Temple Place, Suite 330, Boston, MA 02111-1307 USA.

=== git-rcs-keywords by Martin Turon (turon) at Github ===

MMC includes a pair of git filters (`.git_filters/rcs-keywords.clean` and
`.git_filters/rcs-keywords.smudge`) to automatically update SVN keywords in
`mcx_utils.c`. The two simple filter scripts were licensed under the BSD license
according to this link:

https://github.com/turon/git-rcs-keywords/issues/4

Both filter files were significantly modified by Qianqian Fang.

=== UBJ (UBJSON) library by Steven Braeger ===

* Files: src/ubj/*
* Copyright (c) 2014 Steven Braeger
* URL: https://github.com/Steve132/ubj
* License: MIT License, src/ubj/LICENSE

=== Pybind11 ===

* Files: src/pybind11/*
* Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>
* URL: https://github.com/pybind/pybind11/
* License: BSD-style license, https://github.com/pybind/pybind11/blob/master/LICENSE

== # License ==

MMC is free software, released under the GNU General Public License version 3 or later
(GPLv3+); see `LICENSE.txt`. MMCLAB and PMMC are distributed under the same license.
