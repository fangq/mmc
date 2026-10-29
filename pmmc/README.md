![](https://mcx.space/img/mcx_wiki_banner.png)

# PMMC - Python bindings for Mesh-based Monte Carlo (MMC) photon transport simulator

- Copyright: (C) Qianqian Fang (2025-2026) <q.fang at neu.edu>
- License: GNU Public License V3 or later
- Version: 0.3.12
- URL: https://pypi.org/project/pmmc/
- Github: https://github.com/fangq/mmc

![Linux Python Module](https://github.com/fangq/mmc/actions/workflows/build_linux_manywheel.yml/badge.svg)\
![MacOS Python Module](https://github.com/fangq/mmc/actions/workflows/build_macos_wheel.yml/badge.svg)\
![Windows Python Module](https://github.com/fangq/mmc/actions/workflows/build_windows_wheel.yml/badge.svg)

This module provides a Python binding for Mesh-based Monte Carlo (MMC).
For other binaries, including the standalone executable and the MATLAB bindings,
see [our website](https://mcx.space).

Mesh-based Monte Carlo (MMC) is a 3D Monte Carlo (MC) simulation software for
photon transport in complex turbid media. MMC combines the strengths of the 
MC-based technique and the finite-element (FE) method: on the one hand, it 
can handle general media, including low-scattering ones, as in the MC method; 
on the other hand, it can use an FE-like tetrahedral mesh to represent curved 
boundaries and complex structures, making it even more accurate, flexible, 
and memory efficient. MMC uses the state-of-the-art ray-tracing techniques 
to simulate photon propagation in a mesh space. It has been extensively 
optimized for excellent computational efficiency and portability.


## How to Install

* PIP: ```pip install pmmc```, see https://pypi.org/project/pmmc/

Binary wheels are available for Linux (x86_64), macOS (Apple silicon and Intel) and
Windows. They include the multi-threaded SSE (CPU) and the OpenCL (CPU/GPU) backends.


## Runtime Dependencies
* **CPU or GPU**: An OpenCL-capable CPU or GPU; most modern CPUs or GPUs support OpenCL -
an industrial-standard heterogeneous computing library and specification (https://www.khronos.org/opencl/)
* **OpenCL CPU or GPU runtime/driver**: Both NVIDIA and AMD GPU graphics drivers should contain
out-of-box OpenCL runtimes or drivers; for Intel GPUs, one should install additional OpenCL runtime
support from https://github.com/intel/compute-runtime or install the `intel-opencl-icd` package
if the OS provides (such as Ubuntu 22.04); one can also install an open-source OpenCL runtime
[POCL](http://portablecl.org/), using package manager such as `sudo apt-get install pocl-opencl-icd`. However,
POCL's support is largely limited to CPUs. You **do not need** to install CUDA SDK to use pmmc.
* **Python**: Python 3.6 and newer is required. **Python 2 is not supported**.
* **numpy**: Used to pass/receive volumetric information to/from pmmc. To install, use either conda or pip 
package managers: `pip install numpy` or `conda install numpy`
* **iso2mesh** is an easy-to-use mesh generator for creating the tetrahedral meshed domain
for pmmc, install it with `pip install iso2mesh`
* (optional) **jdata**: Only needed to read/write JNIfTI output files. To install, use pip: `pip install jdata` 
on all operating systems; For Debian-based Linux distributions, you can also install to the system interpreter 
using apt-get: `sudo apt-get install python3-jdata`. See https://pypi.org/project/jdata/ for more details. 
* (optional) **bjdata**: Only needed to read/write BJData/UBJSON files. To install, run `pip install bjdata` 
on all operating systems; For Debian-based Linux distributions, you can also install to the system interpreter 
using apt-get: `sudo apt-get install python3-bjdata`. See https://pypi.org/project/bjdata/ for more details. 
* (optional) **matplotlib**: For plotting the results. To install, run either `pip install matplotlib` or
`conda install matplotlib`

## Build Instructions

### Build Dependencies
* **Operating System**: pmmc and mmc can be compiled on most OSes, including Windows, Linux and MacOS.
* **OpenCL library**: compiling mmc or pmmc requires to link with `libOpenCL.so` on Linux, or `libOpenCL.dylib`
on MacOS or `OpenCL.dll` on Windows. These libraries should have been installed by either graphics driver or
OpenCL runtimes.
* **Python Interpreter**: Python 3.6 or above. The ```pip``` Python package manager and the ```wheel``` package (available
  via ```pip```) are not required but recommended.
* **C/C++ Compiler**: pmmc can be compiled using a wide variety of C compilers, including
  * GNU GCC for Linux, MacOS (installed via MacPorts or brew), and Windows (installed via msys2, mingw64 or cygwin64)
  * Microsoft Visual Studio C/C++ Compiler for Windows.
  * Apple Clang for macOS, available via Xcode.

  Refer to each OS's online documentations for more in-depth information on how to install these compilers.
  MacOS provides built-in OpenCL library support.
* **OpenMP**: The installed C/C++ Compiler should have support for [OpenMP](https://www.openmp.org/). 
  GCC and Microsoft Visual Studio compiler support OpenMP out of the box. Apple Clang, however, requires manual 
  installation of OpenMP libraries for Apple Clang. The easiest way to do this is via the [Brew](https://brew.sh/) package
  manager, preferably after selecting the correct Xcode version:
```bash
brew install libomp
brew link --force libomp
```
* **CMake**: CMake version 3.15 and later is required. Refer to the [CMake website](https://cmake.org/download/) for more information on how to download.
  CMake is also widely available on package managers across all operating systems.

### Build Steps
1. Ensure that ```cmake```, ```python``` and the C/C++ compiler are all located over your ```PATH```.
This can be queried via ```echo $env:PATH``` on Windows or ```echo $PATH``` on Linux. If not, locate them and add their folder to the ```PATH```.

2. Clone the repository and switch to the ```pmmc/``` folder:
```bash
git clone --recursive https://github.com/fangq/mmc.git
cd mmc/pmmc
```
3. One can run `python3 setup.py install` or `python3 -m pip install .` to both locally build and install the module

4. If one only wants to locally build the module, one should run `python3 -m pip wheel .`

5. If the binary module is successfully built locally, you should see a binary wheel file `pmmc-X.X.X-cpXX-cpXX-*.whl`
stored inside the `mmc/pmmc` folder. You can install this wheel package using `python3 -m pip install --force-reinstall pmmc-*.whl`
to force installing this locally compiled `pmmc` module and overwrite any previously installed versions.

6. (optional) To also build the CUDA backend (`cfg["compute"] = "cuda"`, NVIDIA GPUs only), install the CUDA
toolkit and pass `-DBUILD_CUDA=ON` to CMake, for example `CMAKE_ARGS="-DBUILD_CUDA=ON" python3 -m pip install .`


## How to use

The PMMC module is easy to use. You can use the `pmmc.gpuinfo()` function to first list
the OpenCL CPUs and GPUs available on your computer, then call the `run()` function to launch
a photon simulation. The compute backend is selected by `cfg['compute']`: `'opencl'` (default),
`'sse'` (multi-threaded CPU, `cfg['gpuid']=-1`), or `'cuda'` if pmmc was built with CUDA;
`'optix'` is accepted as in mmclab and runs with OpenCL when OptiX is not built in.

A simulation can be defined conveniently in two approaches - a one-liner and a two-liner:

* For the one-liner, one simply passes each MMC simulation setting as a keyword
argument. The supported setting names are compatible with nearly all the input fields
of the MATLAB version of MMC - [MMCLAB](https://github.com/fangq/mmc/blob/master/mmclab/mmclab.m)

```python3
import pmmc
import numpy as np
import matplotlib.pyplot as plt

import iso2mesh as i2m
node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 60], 10, 100)  # create a mesh

gpus = pmmc.gpuinfo()  # list all available GPUs

res = pmmc.run(nphoton=1000000, node=node, elem=elem, elemprop=np.ones(elem.shape[0]), tstart=0, tend=5e-9,
               tstep=5e-9, srcpos=[30,30,0], srcdir=[0,0,1], prop=np.array([[0, 0, 1, 1], [0.005, 1, 0.01, 1.37]]))
res['flux'].shape

plt.imshow(np.log10(res['flux'][30,:, :].squeeze()))
plt.show()
```

* Alternatively, one can also define a Python dict object containing each setting
as a key, and pass on the dict object to `pmmc.run()`

```python3
import pmmc
import numpy as np
cfg = {'nphoton': 1000000, 'node': node, 'elem': elem, 'elemprop': np.ones(elem.shape[0]), 'tstart':0, 'tend':5e-9, 'tstep':5e-9,
       'srcpos': [30,30,0], 'srcdir':[0,0,1], 'prop':[[0,0,1,1],[0.005,1,0.01,1.37]]}
res = pmmc.run(cfg)
```

`pmmc.run()` calls the binary module only. `pmmc.mmclab(cfg)` additionally applies the same
pre- and post-processing as `mmclab.m` (for example, reorienting the elements, finding the
initial element of the source, and unpacking the detected photon and trajectory data into
dicts), and is recommended when porting MMCLAB scripts.
Python ports of all the MMCLAB example scripts can be found in the `example/` folder
(see `example/README.md`).

### Log output

The simulation log (the banner, timings and the absorbed fraction) is controlled by
`cfg['flog']`, as in MATLAB: `2` (default) prints it to Python's `sys.stdout`, so it
shows in Jupyter notebooks; `1` prints it to the stdout of the process, and `0` does the
same without the banner; a string appends it to that file. To silence the log, set
`cfg['flog'] = os.devnull`; the echo of the `cfg` values goes to `sys.stdout` and can
be silenced with `contextlib.redirect_stdout`:

```python3
import os, io, contextlib
cfg['flog'] = os.devnull
with contextlib.redirect_stdout(io.StringIO()):
    res = pmmc.run(cfg)
```

### Reproducible simulations

Photons are scheduled to the CPU threads and GPU threads dynamically for speed, so two runs
with the same `cfg['seed']` agree statistically but do not simulate the identical set of
photons. To reproduce a simulation exactly on the CPU, use the SSE backend with a single
OpenMP thread (set the environment variable `OMP_NUM_THREADS=1` before starting Python).
