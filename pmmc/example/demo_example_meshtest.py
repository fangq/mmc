"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_example_meshtest.m

A sphere in a semi-infinite medium simulated on three meshes of increasing density in one
pmmc.mmclab call (a list of cfg), compared with the analytical solution in
mmc/examples/meshtest/sphdiffsemiinf.mat. This demo requires iso2mesh (pip install iso2mesh).
"""

import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.io import loadmat
from scipy.interpolate import griddata
import pmmc

no, el = i2m.meshasphere([30, 30, 30], 10, 1.0, dovol=False)[:2]  # the surface, as [no, el] = meshasphere()

## shared simulation settings

cfg = {}
cfg["nphoton"] = 3e7
cfg["seed"] = 1648335518
cfg["srcpos"] = [30.0, 30.0, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 1e-10
cfg["prop"] = [[0, 0, 1, 1], [0.002, 1.0, 0.01, 1.37], [0.050, 5.0, 0.9, 1.37]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 0
cfg["method"] = "elem"
cfg = [copy.copy(cfg), copy.copy(cfg), copy.copy(cfg)]

srcpos = [30.0, 30.0, 0.0]


def makemesh(nfull, nodesize, maxvol, boxsize):
    node, elem = i2m.surf2mesh(np.hstack([nfull, np.reshape(nodesize, (-1, 1))]), el, [0, 0, 0], [60.1, 60.1, 60.1], 1, maxvol, regions=np.array([[30, 30, 30]]), dobbx=boxsize)[:2]
    node, elem = i2m.sortmesh(srcpos, node, elem, np.arange(4))[:2]
    elem[:, 4] = elem[:, 4] + 1
    return node, elem


## the three meshes

fixednodes = np.array([[30.0, 30.0, 0.05], [30, 30, 30]])
nfull = np.vstack([no, fixednodes])
cfg[2]["node"], cfg[2]["elem"] = makemesh(nfull, np.concatenate([np.ones(no.shape[0]), [0.5, 3]]), 8, [1.5, 1.5, 1.5, 1.5, 5, 5, 5, 5])
cfg[1]["node"], cfg[1]["elem"] = makemesh(nfull, np.concatenate([np.ones(no.shape[0]), [1, 1]]), 2, [1, 1, 1, 1, 1, 1, 1, 1])
nfull = np.vstack([no, [[30, 30, 30]]])
cfg[0]["node"], cfg[0]["elem"] = makemesh(nfull, np.concatenate([np.ones(no.shape[0]), [3]]), 10, [2, 2, 2, 2, 5, 5, 5, 5])

## run the three simulations

flux = pmmc.mmclab(cfg)

## compare with the analytical solution

xi, yi = np.meshgrid(np.arange(1, 61), np.arange(0, 61))
clines = np.arange(-1, 8.01, 0.5) - 10
ana = loadmat(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../examples/meshtest/sphdiffsemiinf.mat"))
xa, za, phi_ana = ana["xa"], ana["za"], np.array(ana["phi_ana"], dtype=float)
phi_ana[((xa < -12) | (xa > 12)) & (za > -5)] = np.nan
phi_ana[((xa < -10) | (xa > 10)) & (za > 0)] = np.nan
t = np.linspace(0, 2 * np.pi, 200)

fig, axs = plt.subplots(1, 3, figsize=(16, 5))
for k, (ax, plane) in enumerate(zip(axs, (30, 30, 29))):
    s = np.sum(flux[k]["flux"], axis=1) * cfg[0]["tstep"]
    cutpos, cutvalue = i2m.qmeshcut(cfg[k]["elem"][:, :4], cfg[k]["node"], s, np.array([[0, plane, 0], [0, plane, 1], [1, plane, 0]]))[:2]
    cutpos = np.asarray(cutpos)
    vi = griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi))
    ax.contour(xa + 30, za + 31, np.log10(np.abs(phi_ana)), clines, colors=[(0.7, 0.7, 0.7)], linewidths=2)
    ax.contour(np.log10(np.abs(vi)), clines, colors="r", linestyles=":")
    ax.plot(10 * np.cos(t) + 30, 10 * np.sin(t) + 30, "k--", linewidth=2)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title("mesh %d (%d nodes)" % (k, cfg[k]["node"].shape[0]))
plt.show()
