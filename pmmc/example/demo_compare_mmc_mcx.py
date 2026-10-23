"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_compare_mmc_mcx.m

Compare MMC (pmmc) with MCX (pmcx) on a 60x60x60 mm cube with a 20x20x20 mm inclusion, for a
planar source (set planarsrc = 0 for a pencil beam, inclusion = 0 for a homogeneous cube). This
demo requires iso2mesh (pip install iso2mesh) and pmcx (pip install pmcx).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmcx
import pmmc

planarsrc = 1  # set to 1 to compare planar src, set to 0 for pencil beam
inclusion = 1  # set to 1 to add an inclusion in the domain

## prepare the mesh and the MMC input

cfg = {}
cfg["nphoton"] = 3e7
cfg["seed"] = 1648335518

if not inclusion:
    cfg["node"], tmp, cfg["elem"] = i2m.meshabox([0, 0, 0], [60, 60, 60], 2 * np.sqrt(2), 2**3 / 6)
    cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
    cfg["prop"] = [[0, 0, 1, 1], [0.005, 1.0, 0.01, 1.0]]
else:
    node1, face1 = i2m.meshabox([0, 0, 0], [60, 60, 60], 2 * np.sqrt(2), 2**3 / 6)[:2]
    node2, face2 = i2m.meshabox([20, 20, 20], [40, 40, 40], 2 * np.sqrt(2), 2**3 / 6)[:2]
    node, face = i2m.mergemesh(node1, face1[:, :3], node2, face2[:, :3])[:2]
    cfg["node"], cfg["elem"] = i2m.surf2mesh(node, face, None, None, 1, 2**3 / 6, regions=np.array([[1, 1, 1], [30, 30, 30]]), dobbx=0)[:2]
    cfg["elemprop"] = cfg["elem"][:, 4]
    cfg["elem"] = cfg["elem"][:, :4]
    cfg["prop"] = [[0, 0, 1, 1], [0.005, 1.0, 0.01, 1.0], [0.01, 0.2, 0.01, 1.37]]

cfg["srcpos"] = [30.1, 30.2, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 1e-10
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 1
cfg["detpos"] = [[0, 30, 30, 2]]
cfg["method"] = "elem"

if planarsrc:  # when planar src is used, re-tessellate the mesh
    cfg["srctype"] = "planar"
    cfg["srcpos"] = [10, 10, -10]
    cfg["srcparam1"] = [40, 0, 0, 0]
    cfg["srcparam2"] = [0, 40, 0, 0]
    srcdef = {k: cfg[k] for k in ("srctype", "srcpos", "srcdir", "srcparam1", "srcparam2")}
    sdom = pmmc.mmcsrcdomain(srcdef, [cfg["node"].min(axis=0), cfg["node"].max(axis=0)])[0]
    cfg["node"], cfg["elem"] = pmmc.mmcaddsrc(cfg["node"], np.hstack([cfg["elem"][:, :4], np.reshape(cfg["elemprop"], (-1, 1))]), sdom)
    cfg["elemprop"] = cfg["elem"][:, 4]
    cfg["elem"] = cfg["elem"][:, :4]

## run MMC

cfg = pmmc.mmclab(cfg, "prep")
cube = pmmc.mmclab(cfg)["flux"]

## run MCX on the same domain

cfgx = {k: v for k, v in cfg.items() if k not in ("node", "elem", "elemprop", "debuglevel", "facenb", "evol", "e0", "isreoriented", "method")}
dim = 60
cfgx["vol"] = np.ones([dim, dim, dim], dtype="uint8")
cfgx["faststep"] = 0
if inclusion:
    cfgx["vol"][20:40, 20:40, 20:40] = 2
cfgx["nphoton"] = cfgx["vol"].size / cfg["node"].shape[0] * cfg["nphoton"]  # for mcx, nphoton scaled by the voxel/node ratio
cfgx["autopilot"] = 1
cfgx["gpuid"] = 1
cfgx["issrcfrom0"] = 1  # treating the lower-outer corner of the grid as [0 0 0]
cfgx["issavedet"] = 0  # mcxlab with one output
cubex = pmcx.mcxlab(cfgx)["flux"]

## compare the CW fluence at y=30.5

cwcb = np.sum(cube, axis=1)
cwcbx = np.sum(cubex, axis=3)

cutpos, cutvalue = i2m.qmeshcut(cfg["elem"][:, :4], cfg["node"], cwcb, "y=30.5")[:2]
cutpos = np.asarray(cutpos)
xi, yi = np.meshgrid(np.arange(0.5, 60), np.arange(0.5, 60))  # interpolate mmc solution at voxel centers
vi = griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi))
vix = cwcbx[:, 29, :].T

clines = np.arange(-8, -1.49, 0.25)
plt.figure()
if not planarsrc and not inclusion:
    dets = np.column_stack([xi.ravel(), 30.2 * np.ones(xi.size), yi.ravel()])
    phicw = np.reshape(pmmc.cwdiffusion(0.005, 1.0, 0, np.array([[30.1, 30.2, 0]]), dets)[0], xi.shape)
    plt.contour(xi, yi, np.log10(np.maximum(phicw, 1e-8)), np.arange(-4, -1.49, 0.5), colors="k")
plt.contour(xi, yi, np.log10(np.abs(vi * cfg["tstep"])), clines, colors="r", linestyles=":")
plt.contour(xi, yi, np.log10(np.abs(vix * cfg["tstep"])), clines, colors="b", linestyles="--")
if inclusion:
    plt.plot([20, 20, 40, 40, 20], [20, 40, 40, 20, 20], "k--")
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.legend(["MMC", "MCX"] if planarsrc or inclusion else ["Diffusion", "MMC", "MCX"])
plt.title("CW fluence at y=30.5, red: MMC, blue: MCX")
plt.show()
