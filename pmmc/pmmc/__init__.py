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

"""PMMC - Python bindings for Monte Carlo eXtreme (OpenCL) photon transport simulator

Example usage:

# To list available GPUs
import pmmc
pmmc.gpuinfo()

# To generate tetrahedral mesh
import pyvista as pv
import tetgen
import numpy as np

box = pv.Box(bounds=(0, 60, 0, 60, 0, 60))
box_tri = box.triangulate()
tet = tetgen.TetGen(box_tri)
node, elem = tet.tetrahedralize(order=1, minratio=1.5, mindihedral=20, switches='pq1.2a50')
elem = elem + 1

# To run a simulation
cfg = {'nphoton': 1000000, 'node': node, 'elem': elem, 'elemprop': np.ones(elem.shape[0]),
       'tstart':0, 'tend':5e-9, 'tstep':5e-9, 'srcpos': [30,30,0], 'srcdir':[0,0,1],
       'prop':[[0,0,1,1],[0.005,1,0.01,1.37]]}
res = pmmc.run(cfg)
"""

try:
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
    )
except ImportError:  # pragma: no cover
    print(
        "please first install pmcx module to use utility functions such as 'detweight', 'meanpath' etc"
    )

__version__ = "0.2.0"

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
)
