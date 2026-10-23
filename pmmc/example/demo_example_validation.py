"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_example_validation.m

Validate MMC against the analytical diffusion solution of a semi-infinite medium: a 60x60x60 mm
cube (5 tetrahedra per 2 mm cube) with a pencil beam on the top surface. This demo requires
iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmmc

## prepare simulation input

cfg = {}
cfg["nphoton"] = 3e6
cfg["seed"] = 1648335518
cfg["node"], cfg["elem"] = i2m.meshgrid5(np.arange(0, 61, 2), np.arange(0, 61, 2), np.arange(0, 61, 2))[:2]  # genT5mesh in MATLAB
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [30.1, 30.2, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 1e-10
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1.0, 0.01, 1.0]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 0
cfg["method"] = "elem"

## run the simulation

cube = pmmc.mmclab(cfg)["flux"]

## plot the time-domain response at one node and the CW fluence on y=30.2

c0 = 299792458000
twin = np.arange(cfg["tstart"] + cfg["tstep"] / 2, cfg["tend"], cfg["tstep"])
gates = len(twin)

cwcb = np.sum(cube, axis=1)
cutpos, cutvalue = i2m.qmeshcut(cfg["elem"][:, :4], cfg["node"], cwcb, np.array([[0, 30.2, 0], [0, 30.2, 1], [1, 30.2, 0]]))[:2]
cutpos = np.asarray(cutpos)
xi, yi = np.meshgrid(np.arange(0, 61, 2), np.arange(0, 61, 2))
vi = griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi))

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

srcpos = [30.1, 30.2, 0]
detpos = [30, 14, 10]
axs[0].semilogy(np.arange(1, gates + 1) / 10, np.ravel(pmmc.tddiffusion(0.005, 1, c0, 0, srcpos, detpos, twin)), "r")
idx = np.where(np.all(np.isclose(cfg["node"], detpos), axis=1))[0]
axs[0].semilogy(np.arange(1, gates + 1) / 10, cube[idx[0], :], "+")
axs[0].set_xlabel("t (ns)")
axs[0].set_ylabel("fluence rate (1/(mm^2 s))")
axs[0].legend(["Diffusion", "MMC"])
axs[0].set_title("time-domain fluence rate at [30, 14, 10]")

clines = np.arange(-1.5, -8.01, -0.5)[::-1]
dets = np.column_stack([xi.ravel(), 30.2 * np.ones(xi.size), yi.ravel()])
phicw = np.reshape(pmmc.cwdiffusion(0.005, 1.0, 0, np.array([srcpos]), dets)[0], xi.shape)  # [Phi, r] in MATLAB
axs[1].contour(xi, yi, np.log10(np.maximum(phicw, 1e-8)), np.arange(-4, -1.49, 0.5), colors="k")
axs[1].contour(xi, yi, np.log10(np.abs(vi * cfg["tstep"])), clines, colors="r", linestyles=":")
axs[1].invert_yaxis()
axs[1].set_aspect("equal")
axs[1].set_title("CW fluence at y=30.2 (black: diffusion, red: MMC)")
plt.show()
