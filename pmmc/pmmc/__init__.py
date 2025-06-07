# Copyright (c) 2022-2023 Matin Raayai Ardakani <raayaiardakani.m at northeastern.edu>. All rights reserved.
# Copyright (c) 2022-2024 Qianqian Fang <q.fang at neu.edu>. All rights reserved.
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

box = pv.Box(bounds=(0, 60, 0, 60, 0, 60))
box_tri = box.triangulate()
tet = tetgen.TetGen(box_tri)
tet.tetrahedralize(order=1, minratio=1.5, mindihedral=20)
tetra_mesh = tet.grid
node = tetra_mesh.points
elem = tetra_mesh.cells_dict[pv.CellType.TETRA] + 1

# To run a simulation
res = pmmc.run(nphoton=1000000, node=node, elem=elem, elemprop=np.ones(elem.shape[0]),
               tstart=0, tend=5e-9, tstep=5e-9, srcpos=[30,30,0],
               srcdir=[0,0,1], prop=np.array([[0, 0, 1, 1], [0.005, 1, 0.01, 1.37]]))
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

__version__ = "0.1.1"

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
