"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmcl_b1_b1d.m

Benchmarks B1 and B1D of the MMC paper: a homogeneous 60 mm cube simulated with regular MMC
(method='elem', 5 tetrahedra per 2 mm cube) and with dual-grid MMC (method='grid', a 6-element
mesh of the cube), each on the GPU (gpuid=1) and on the CPU (gpuid=-1). This demo requires
iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmmc

## B1: regular MMC on a fine mesh

cfg = {}
cfg["nphoton"] = 1e6
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

cfg["gpuid"] = 1
b1gpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=1)
cfg["gpuid"] = -1
b1cpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=1)

xi, yi = np.meshgrid(np.arange(0.5, 60), np.arange(0.5, 60))
plane = np.array([[0, 30.5, 0], [0, 30.5, 1], [1, 30.5, 0]])
phi = []
for b1 in (b1cpu, b1gpu):
    cutpos, cutvalue = i2m.qmeshcut(cfg["elem"][:, :4], cfg["node"], b1, plane)[:2]
    cutpos = np.asarray(cutpos)
    phi.append(griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi)))
Phicpu, Phigpu = phi

## B1D: dual-grid MMC on a 6-element mesh of the cube

cfg["node"], cfg["elem"] = i2m.meshgrid6([0, 60], [0, 60], [0, 60])[:2]
cfg["elem"][:, :4] = i2m.meshreorient(cfg["node"], cfg["elem"][:, :4].copy())[0]
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["method"] = "grid"

cfg["gpuid"] = 1
Phidgpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=3)
cfg["gpuid"] = -1
Phidcpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=3)

## compare the four CW fluence maps on y=30.5

clines = np.arange(-8, 0.01, 0.5)
plt.figure()
plt.contourf(xi, yi, np.log10(np.maximum(Phicpu * cfg["tstep"], 1e-8)), clines)
plt.contour(xi, yi, np.log10(np.abs(Phigpu * cfg["tstep"])), clines, colors="g", linestyles="--", linewidths=2)
plt.contour(xi, yi, np.log10(np.abs(Phidcpu[:-1, 30, :-1].T * cfg["tstep"])), clines, colors="b", linestyles="--", linewidths=2)
plt.contour(xi, yi, np.log10(np.abs(Phidgpu[:-1, 30, :-1].T * cfg["tstep"])), clines, colors="r", linewidths=1)
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("B1 CPU (filled), B1 GPU (green), B1D CPU (blue), B1D GPU (red)")
plt.show()
