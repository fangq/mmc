"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_dmmc_sphshells.m

Dual-grid MMC (DMMC, method='grid') vs. regular MMC (method='elem') vs. MCX on nested spherical
shells (a CSF-like thin layer between two tissues). This demo requires iso2mesh
(pip install iso2mesh) and pmcx (pip install pmcx).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmcx
import pmmc


def sphere(r, tsize):
    node, face = i2m.meshasphere([30, 30, 30], r, tsize, dovol=False)[:2]  # the surface, as [node, face] = meshasphere()
    return i2m.removeisolatednode(node, face)[:2]


nsph_10, fsph_10 = sphere(10, 1.0)
nsph_23, fsph_23 = sphere(23, 2.0)
nsph_25, fsph_25 = sphere(25, 2.0)

## shared simulation settings

cfg = {}
cfg["nphoton"] = 3e6
cfg["seed"] = 1648335518
cfg["srcpos"] = [30.5, 30.5, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["prop"] = [[0, 0, 1, 1], [0.02, 7.0, 0.89, 1.37], [0.004, 0.009, 0.89, 1.37], [0.02, 9.0, 0.89, 1.37], [0.05, 0.0, 1.0, 1.37]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 1
cfgs = [copy.copy(cfg), copy.copy(cfg)]

## DMMC: a coarse mesh of the spheres inside a box

nbox, ebox = i2m.meshgrid6([0, 60], [0, 60], [0, 60])[:2]
fbox = i2m.volface(ebox)[0]
no, fc = i2m.mergemesh(nsph_10, fsph_10, nsph_23, fsph_23, nsph_25, fsph_25, nbox, fbox)[:2]
no, fc = i2m.removeisolatednode(no, fc[:, :3])[:2]
cfgs[0]["node"], cfgs[0]["elem"] = i2m.surf2mesh(
    no, fc, [0, 0, 0], [60.1, 60.1, 60.1], 1, 100, regions=np.array([[1, 1, 1], [30, 30, 6], [30, 30, 15], [30, 30, 30]]), cmdopt="-A -q -Y"
)[:2]  # thin layer
cfgs[0]["method"] = "grid"

## MMC: a dense mesh of the spheres

no, fc = i2m.mergemesh(nsph_10, fsph_10, nsph_23, fsph_23, nsph_25, fsph_25)[:2]
no, fc = i2m.removeisolatednode(no, fc[:, :3])[:2]
srcpos = [30.5, 30.5, 0.0]
cfgs[1]["node"], cfgs[1]["elem"] = i2m.surf2mesh(
    no, fc, [0, 0, 0], [60, 60, 60], 1, 0.25, regions=np.array([[1, 1, 1, 0.1], [30, 30, 6, 0.1], [30, 30, 17, 0.1], [30, 30, 30, 0.1]]), dobbx=[1, 1, 1, 1, 1, 1, 1, 1]
)[:2]
cfgs[1]["node"], cfgs[1]["elem"] = i2m.sortmesh(srcpos, cfgs[1]["node"], cfgs[1]["elem"], np.arange(4))[:2]
cfgs[1]["method"] = "elem"

phimmc = pmmc.mmclab(cfgs)

## MCX on the voxelized domain

cfg = {}
cfg["seed"] = 0x623F9A9E
cfg["nphoton"] = 1e8
dim = 60
xi, yi, zi = np.meshgrid(*([np.arange(0.5, dim)] * 3), indexing="ij")
dist = (xi - 30) ** 2 + (yi - 30) ** 2 + (zi - 30) ** 2
vol = np.ones(xi.shape)
vol[dist < 625] = 2
vol[dist < 529] = 3
vol[dist < 100] = 4
cfg["vol"] = vol.astype("uint8")
cfg["srcpos"] = [30.5, 30.5, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["issrcfrom0"] = 1
cfg["prop"] = [[0, 0, 1, 1], [0.02, 7.0, 0.89, 1.37], [0.004, 0.009, 0.89, 1.37], [0.02, 9.0, 0.89, 1.37], [0.05, 0.0, 1.0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["autopilot"] = 1
cfg["gpuid"] = 1
cfg["isreflect"] = 1  # disable reflection at exterior boundaries
cfg["unitinmm"] = 1  # resolution
cfg["issavedet"] = 0  # mcxlab with one output
print("running MCX simulation ...")
phimcx = pmcx.mcxlab(cfg)
phi_mcx = np.sum(phimcx["flux"], axis=3)

## compare on y=30.5

xx, yy = np.meshgrid(np.arange(0.5, 60), np.arange(0.5, 60))
s1 = np.sum(phimmc[1]["flux"], axis=1)
cutpos, cutvalue = i2m.qmeshcut(cfgs[1]["elem"][:, :4], cfgs[1]["node"], s1, np.array([[0, 30.5, 0], [0, 30.5, 1], [1, 30.5, 0]]))[:2]
cutpos = np.asarray(cutpos)
phi_mmc = griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xx, yy))
phi_dmmc = np.sum(phimmc[0]["flux"], axis=3)[:60, :60, :60]

clines = np.arange(-10, 10.01, 0.5)
plt.figure()
plt.contourf(np.log10(np.abs(phi_mmc)), clines)
plt.contour(np.log10(np.abs(phi_dmmc[:, 30, :].T)), clines, linestyles="--", colors="w", linewidths=1.5)
plt.contour(np.log10(np.abs(phi_mcx[:, 30, :].T)), clines, linestyles="--", colors=[(0.91, 0.41, 0.17)], linewidths=1.5)
t = np.linspace(0, 2 * np.pi, 200)
for r in (10, 23, 25):
    plt.plot(r * np.cos(t) + 30.5, r * np.sin(t) + 30.5, "k--", linewidth=2)
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("CW fluence, filled: MMC, white: DMMC, orange: MCX")
plt.show()
