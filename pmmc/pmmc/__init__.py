# Copyright (c) 2025 Qianqian Fang <q.fang at neu.edu>. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""PMMC - Python bindings for Mesh-based Monte Carlo (MMC) photon transport simulator

Example usage:

# To list available GPUs
import pmmc
pmmc.gpuinfo()

import numpy as np

# To generate tetrahedral mesh
import iso2mesh as i2m
node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 60], 10, 100)  # create a mesh

# To run a simulation
cfg = {'nphoton': 1000000, 'node': node, 'elem': elem, 'elemprop': np.ones(elem.shape[0]),
       'tstart':0, 'tend':5e-9, 'tstep':5e-9, 'srcpos': [30,30,0], 'srcdir':[0,0,1],
       'prop':[[0,0,1,1],[0.005,1,0.01,1.37]]}
res = pmmc.run(cfg)
"""

import sys

try:
    if sys.platform.startswith("win"):
        import os, sparse_numba, ctypes

        ctypes.CDLL(
            os.path.join(
                os.path.dirname(sparse_numba.__file__),
                "vendor",
                "superlu",
                "bin",
                "libgomp-1.dll",
            )
        )

    from _pmmc import gpuinfo, run, version
except ImportError:  # pragma: no cover
    print("the pmmc binary extension (_pmmc) is not compiled! please compile first")

try:
    from pmcx import (
        detweight,
        cwdref,
        meanpath,
        meanscat,
        dettpsf,
        dettime,
        tddiffusion,
        getdistance,
        detphoton,
        mcxlab,
        cwdiffusion,
        cwfluxdiffusion,
        cwfluencediffusion,
        dcsg1,
        mcxcreate,
        rfreplay,
        rfmusreplay,
        loadmc2,
        loadmch,
        loadfile,
        mcx2json,
        json2mcx,
        loadnii,
        preview,
        plotshapes,
        plotphotons,
        plotvol,
    )
except ImportError:  # pragma: no cover
    print(
        "please first install pmcx module to use utility functions such as 'detweight', 'meanpath' etc"
    )

__version__ = "0.3.5"

__all__ = (
    "gpuinfo",
    "run",
    "version",
    "detweight",
    "cwdref",
    "meanpath",
    "meanscat",
    "dettpsf",
    "dettime",
    "tddiffusion",
    "getdistance",
    "detphoton",
    "mcxlab",
    "cwdiffusion",
    "cwfluxdiffusion",
    "cwfluencediffusion",
    "dcsg1",
    "mcxcreate",
    "rfreplay",
    "rfmusreplay",
    "loadmc2",
    "loadmch",
    "loadfile",
    "mcx2json",
    "json2mcx",
    "loadnii",
    "preview",
    "plotshapes",
    "plotphotons",
    "plotvol",
)
