"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmcl_b4.m

Benchmark B4 of the MMC paper: the mcxyz skin-vessel model (a blood vessel in a 3-layer tissue)
with a disk source, simulated on the GPU (gpuid=1) and on the CPU (gpuid=-1). This demo requires
iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmmc

## create the mesh: a vessel merged with a 3-layer slab

no, fc = i2m.latticegrid([0, 200], [0, 200], [20, 32, 200])[:2]  # create a 3-layer tissue
fc2 = np.array(fc)
fc = np.vstack([fc2[:, [0, 1, 2]], fc2[:, [0, 2, 3]]])
ncy, fcy = i2m.meshacylinder([-1, 99.5, 99.5], [201, 99.5, 99.5], 20, 5)[:2]  # add the vessel
newnode, newelem = i2m.surfboolean(no, fc, "first", ncy, fcy[:, :3])[:2]  # merge the two domains

c0 = np.array([10, 100, 150, 26])
seeds = np.column_stack([np.ones((4, 2)) * 100, c0])  # define the regions by index

cfg_mmc = {}
cfg_mmc["node"], cfg_mmc["elem"] = i2m.s2m(newnode, newelem[:, :3], 1, 4, "tetgen", seeds, None)[:2]  # creating the merged mesh domain
voxellen = 0.005
cfg_mmc["node"] = cfg_mmc["node"] * voxellen
cfg_mmc["method"] = "elem"
cfg_mmc["elemprop"] = cfg_mmc["elem"][:, 4]
cfg_mmc["elem"] = cfg_mmc["elem"][:, :4]
cfg_mmc["prop"] = [[0.0000, 0.0, 1.0000, 1], [3.5640e-05, 1.0000, 1.0000, 1.3700], [23.0543, 9.3985, 0.9000, 1.3700], [0.0458, 35.6541, 0.9000, 1.3700], [1.6572, 37.5940, 0.9000, 1.3700]]

## define the disk source and add its domain to the mesh

cfg_mmc["srcpos"] = list(np.array([100, 100, -3]) * voxellen)
cfg_mmc["srcdir"] = [0, 0, 1]
cfg_mmc["tstart"] = 0
cfg_mmc["tend"] = 5e-8
cfg_mmc["tstep"] = 5e-8
cfg_mmc["outputtype"] = "flux"
cfg_mmc["minenergy"] = 0.01
cfg_mmc["srctype"] = "disk"
cfg_mmc["srcparam1"] = [0.3, 0, 0, 0]

srcdef = {"srctype": cfg_mmc["srctype"], "srcpos": cfg_mmc["srcpos"], "srcdir": cfg_mmc["srcdir"], "srcparam1": cfg_mmc["srcparam1"], "srcparam2": []}
sdom = pmmc.mmcsrcdomain(srcdef, [cfg_mmc["node"].min(axis=0), cfg_mmc["node"].max(axis=0)])[0]
cfg_mmc["node"], cfg_mmc["elem"] = pmmc.mmcaddsrc(cfg_mmc["node"], np.hstack([cfg_mmc["elem"], cfg_mmc["elemprop"][:, None]]), sdom)
cfg_mmc["elemprop"] = cfg_mmc["elem"][:, 4]
cfg_mmc["elem"] = cfg_mmc["elem"][:, :4]

## run on the GPU and on the CPU

cfg_mmc["nphoton"] = 1e8
cfg_mmc["seed"] = 1648335518
cfg_mmc["debuglevel"] = "TP"
cfg_mmc["isreflect"] = 0

cfg_mmc["gpuid"] = 1
b4gpu = pmmc.mmclab(cfg_mmc)["flux"] * cfg_mmc["tstep"]
cfg_mmc["gpuid"] = -1
b4cpu = pmmc.mmclab(cfg_mmc)["flux"] * cfg_mmc["tstep"]

## compare on the x=0.4975 mm plane

yi, zi = np.meshgrid(np.arange(0.0025, 0.998, 0.005), np.arange(0.1025, 0.998, 0.005), indexing="ij")
plane = np.array([[0.4975, 0, 0], [0.4975, 0, 1], [0.4975, 1, 0]])
interp = []
for b4 in (b4cpu, b4gpu):
    cutpos, cutvalue = i2m.qmeshcut(cfg_mmc["elem"][:, :4], cfg_mmc["node"], np.ravel(b4[:, 0]), plane)[:2]
    cutpos = np.asarray(cutpos)
    interp.append(griddata((cutpos[:, 1], cutpos[:, 2]), np.ravel(cutvalue), (yi, zi)))

yyy, zzz = np.meshgrid(np.arange(-0.4975, 0.498, 0.005), np.arange(0.1025, 0.998, 0.005), indexing="ij")
clines = np.arange(-2, 1.01, 0.25)
plt.figure()
plt.contourf(yyy, zzz, np.log10(interp[0]), clines)
plt.contour(yyy, zzz, np.log10(interp[1]), clines, colors=[(0.91, 0.41, 0.17)], linestyles=":", linewidths=2)
t = np.linspace(0, 2 * np.pi, 11)
plt.plot(0.1 * np.cos(t) - 0.0025, 0.1 * np.sin(t) + 0.4975, "k--", linewidth=2)
plt.plot([-0.5, 0.5], [0.16, 0.16], "k--", linewidth=2)
plt.xlim([-0.5, 0.5])
plt.ylim([0.1, 1])
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("B4 CW fluence (log10), filled: CPU, orange: GPU")
plt.show()
