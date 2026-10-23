"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mcxyz_skinvessel.m

The skin-vessel benchmark of mcxyz (https://omlc.org/software/mc/mcxyz/): a blood vessel in a
3-layer tissue (air, epidermis, dermis), illuminated by a disk source. This demo requires
iso2mesh (pip install iso2mesh).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshplot import plotcut

## create the mesh: a vessel merged with a 3-layer slab

no, fc = i2m.latticegrid([0, 200], [0, 200], [20, 32, 200])[:2]  # create a 3-layer tissue
no[-1, :] = no[-1, :] + 1e-5
fc2 = np.array(fc)
fc = np.vstack([fc2[:, [0, 1, 2]], fc2[:, [0, 2, 3]]])
ncy, fcy = i2m.meshacylinder([-1, 99.5, 99.5], [201, 99.5, 99.5], 20, 5)[:2]  # add the vessel
newnode, newelem = i2m.surfboolean(no, fc, "first", ncy, fcy[:, :3])[:2]  # merge the two domains

c0 = np.array([10, 100, 150, 26])
seeds = np.column_stack([np.ones((4, 2)) * 100, c0])  # define the regions by index

cfg = {}
cfg["node"], cfg["elem"] = i2m.s2m(newnode, newelem[:, :3], 1, 30, "tetgen", seeds, None)[:2]  # creating the merged mesh domain
cfg["unitinmm"] = 0.005
cfg["method"] = "elem"
cfg["elemprop"] = cfg["elem"][:, 4]
cfg["elem"] = cfg["elem"][:, :4]
cfg["prop"] = [[0.0000, 0.0, 1.0000, 1], [3.5640e-05, 1.0000, 1.0000, 1.3700], [23.0543, 9.3985, 0.9000, 1.3700], [0.0458, 35.6541, 0.9000, 1.3700], [1.6572, 37.5940, 0.9000, 1.3700]]

## define the disk source and add its domain to the mesh

cfg["srcpos"] = [100, 100, -1]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-8
cfg["tstep"] = 5e-8
cfg["outputtype"] = "flux"
cfg["minenergy"] = 0.01
cfg["srctype"] = "disk"
cfg["srcparam1"] = list(np.array([0.3, 0, 0, 0]) / cfg["unitinmm"])  # in grid unit

srcdef = {"srctype": cfg["srctype"], "srcpos": cfg["srcpos"], "srcdir": cfg["srcdir"], "srcparam1": cfg["srcparam1"], "srcparam2": []}
sdom = pmmc.mmcsrcdomain(srcdef, [cfg["node"].min(axis=0), cfg["node"].max(axis=0)])[0]
cfg["node"], cfg["elem"] = pmmc.mmcaddsrc(cfg["node"], np.hstack([cfg["elem"], cfg["elemprop"][:, None]]), sdom)
cfg["elemprop"] = cfg["elem"][:, 4]
cfg["elem"] = cfg["elem"][:, :4]

## run the simulation

cfg["nphoton"] = 1e7
cfg["seed"] = 1648335518
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 0

flux = pmmc.mmclab(cfg)["flux"]
fluxcw = np.sum(flux, axis=1) * cfg["tstep"] * 100  # multiplying 100 converts 1/mm^2 to 1/cm^2 to match mcxyz output

## plot the fluence in the tissue on the x=0.5 mm plane

plt.figure()
plotcut(cfg["node"] * cfg["unitinmm"], cfg["elem"][cfg["elemprop"] > 0], np.log10(fluxcw), "x=0.5", title="log10 CW fluence (1/cm^2) at x=0.5 mm")
plt.gca().invert_yaxis()
plt.show()
