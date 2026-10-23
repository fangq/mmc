"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_dualmesh_output.m

Dual-grid MMC (DMMC, method='grid', a coarse mesh with a voxel output grid) vs. regular MMC
(method='elem', a dense mesh) vs. MCX for a sphere in a box. This demo requires iso2mesh
(pip install iso2mesh) and pmcx (pip install pmcx).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmcx
import pmmc

nsph, fsph = i2m.meshasphere([30, 30, 30], 10, 1.0, dovol=False)[:2]  # the surface, as [node, face] = meshasphere()
nsph, fsph = i2m.removeisolatednode(nsph, fsph)[:2]

## shared simulation settings

cfg = {}
cfg["nphoton"] = 1e7
cfg["seed"] = 1648335518
cfg["srcpos"] = [30, 30.001, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["prop"] = [[0, 0, 1, 1], [0.002, 1.0, 0.01, 1.37], [0.050, 0, 1, 1.37]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 1
cfgs = [copy.copy(cfg), copy.copy(cfg)]

## DMMC: a coarse mesh of the sphere in a box

nbox, ebox = i2m.meshgrid6([0, 60], [0, 60], [0, 60])[:2]
fbox = i2m.volface(ebox)[0]
no, fc = i2m.mergemesh(nbox, fbox, nsph, fsph)[:2]
cfgs[0]["node"], cfgs[0]["elem"] = i2m.surf2mesh(no, fc[:, :3], [0, 0, 0], [60.1, 60.1, 60.1], 1, 100, regions=np.array([[1, 1, 1], [30, 30, 30]]), cmdopt="-Y -A")[:2]
cfgs[0]["method"] = "grid"

## MMC: a dense mesh with fixed nodes near the source

no, el = i2m.meshasphere([30, 30, 30], 10, 1.0, dovol=False)[:2]
nodesize = np.concatenate([np.ones(no.shape[0]), [1, 1]])
srcpos = [30.0, 30.001, 0.0]
fixednodes = np.array([[30.0, 30.0, 0.05], [30, 30, 30]])
nfull = np.vstack([no, fixednodes])
cfgs[1]["node"], cfgs[1]["elem"] = i2m.surf2mesh(np.hstack([nfull, nodesize[:, None]]), el, [0, 0, 0], [60.1, 60.1, 60.1], 1, 2, regions=np.array([[30, 30, 30]]), dobbx=[1, 1, 1, 1, 1, 1, 1, 1])[:2]
cfgs[1]["node"], cfgs[1]["elem"] = i2m.sortmesh(srcpos, cfgs[1]["node"], cfgs[1]["elem"], np.arange(4))[:2]
cfgs[1]["elem"][:, 4] = cfgs[1]["elem"][:, 4] + 1
cfgs[1]["method"] = "elem"

phimmc = pmmc.mmclab(cfgs)

## MCX on the voxelized domain

cfg = {}
cfg["seed"] = 0x623F9A9E
cfg["nphoton"] = 3e7
dim = 60
xi, yi, zi = np.meshgrid(*([np.arange(0.5, dim)] * 3), indexing="ij")
dist = (xi - 30) ** 2 + (yi - 30) ** 2 + (zi - 30) ** 2
vol = np.ones(xi.shape)
vol[dist < 100] = 2
cfg["vol"] = vol.astype("uint8")
cfg["srcpos"] = [30, 30.001, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["issrcfrom0"] = 1
cfg["prop"] = [[0, 0, 1, 1], [0.002, 1.0, 0.01, 1.37], [0.050, 0, 1, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["autopilot"] = 1
cfg["gpuid"] = 1
cfg["isreflect"] = 1
cfg["issavedet"] = 0  # mcxlab with one output
print("running MCX simulation ...")
phimcx = pmcx.mcxlab(cfg)

## compare on y=30.5

clines = np.arange(-7, -0.99, 0.25)
xx, zz = np.meshgrid(np.arange(0.5, 60), np.arange(0.5, 60))
CWfluence_mmc = np.sum(phimmc[1]["flux"], axis=1) * cfg["tstep"]
cutpos, cutvalue = i2m.qmeshcut(cfgs[1]["elem"][:, :4], cfgs[1]["node"], CWfluence_mmc, np.array([[0, 30.5, 0], [0, 30.5, 1], [1, 30.5, 0]]))[:2]
cutpos = np.asarray(cutpos)
CWfluence_mmc_grid = griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xx, zz))
CWfluence_dmmc = np.sum(phimmc[0]["flux"], axis=3) * cfg["tstep"]
CWfluence_mcx = np.sum(phimcx["flux"], axis=3) * cfg["tstep"]

plt.figure()
plt.contourf(xx, zz, np.log10(np.abs(CWfluence_mmc_grid)), clines)
plt.contour(xx, zz, np.log10(np.abs(CWfluence_dmmc[:-1, 30, :-1].T)), clines, linestyles="--", colors="w", linewidths=1.5)
plt.contour(xx, zz, np.log10(np.abs(CWfluence_mcx[:, 30, :].T)), clines, linestyles="--", colors=[(0.91, 0.41, 0.17)], linewidths=1.5)
t = np.linspace(0, 2 * np.pi, 200)
plt.plot(10 * np.cos(t) + 30, 10 * np.sin(t) + 30, "k--", linewidth=2)
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("CW fluence, filled: MMC, white: DMMC, orange: MCX")
plt.show()
