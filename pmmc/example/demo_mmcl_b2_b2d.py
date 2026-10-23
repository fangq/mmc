"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmcl_b2_b2d.m

Benchmarks B2 and B2D of the MMC paper: nested spherical shells with a thin CSF-like layer,
simulated with regular MMC (method='elem') and dual-grid MMC (method='grid'), each on the GPU and
the CPU, and compared with MCX-CL (pmcxcl). This demo requires iso2mesh (pip install iso2mesh)
and pmcxcl (pip install pmcxcl).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmcxcl
import pmmc


def sphere(r, tsize):
    node, face = i2m.meshasphere([30, 30, 30], r, tsize, dovol=False)[:2]  # the surface, as [node, face] = meshasphere()
    return i2m.removeisolatednode(node, face)[:2]


nsph_10, fsph_10 = sphere(10, 1.0)
nsph_23, fsph_23 = sphere(23, 2.0)
nsph_25, fsph_25 = sphere(25, 2.0)

## shared simulation settings

cfg = {}
cfg["nphoton"] = 1e8
cfg["seed"] = 1648335518
cfg["srcpos"] = [30, 30.1, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["prop"] = [[0, 0, 1, 1], [0.02, 7.0, 0.89, 1.37], [0.0004, 0.009, 0.89, 1.37], [0.02, 9.0, 0.89, 1.37], [0.05, 0.0, 1.0, 1.37]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 1
cfgs = [copy.copy(cfg), copy.copy(cfg)]

## B2D mesh: the spheres in a box, for dual-grid MMC

nbox, ebox = i2m.meshgrid6([0, 60], [0, 60], [0, 60])[:2]
fbox = i2m.volface(ebox)[0]
no, fc = i2m.mergemesh(nsph_10, fsph_10, nsph_23, fsph_23, nsph_25, fsph_25, nbox, fbox)[:2]
no, fc = i2m.removeisolatednode(no, fc[:, :3])[:2]
cfgs[0]["node"], cfgs[0]["elem"] = i2m.surf2mesh(
    no, fc, [0, 0, 0], [60.1, 60.1, 60.1], 1, 100, regions=np.array([[1, 1, 1], [30, 30, 6], [30, 30, 15], [30, 30, 30]]), cmdopt="-A -q -Y"
)[:2]  # thin layer
cfgs[0]["method"] = "grid"

## B2 mesh: a dense mesh of the spheres

no, fc = i2m.mergemesh(nsph_10, fsph_10, nsph_23, fsph_23, nsph_25, fsph_25)[:2]
no, fc = i2m.removeisolatednode(no, fc[:, :3])[:2]
srcpos = cfgs[0]["srcpos"]
cfgs[1]["node"], cfgs[1]["elem"] = i2m.surf2mesh(
    no, fc, [0, 0, 0], [60, 60, 60], 1, 0.25, regions=np.array([[1, 1, 1, 0.1], [30, 30, 6, 0.1], [30, 30, 17, 0.1], [30, 30, 30, 0.1]]), dobbx=[1, 1, 1, 1, 1, 1, 1, 1]
)[:2]
cfgs[1]["node"], cfgs[1]["elem"] = i2m.sortmesh(srcpos, cfgs[1]["node"], cfgs[1]["elem"], np.arange(4))[:2]
cfgs[1]["method"] = "elem"

## run B2D and B2 on the GPU and the CPU

cfgs[0]["gpuid"] = 1
b2dgpu = np.sum(pmmc.mmclab(cfgs[0])["flux"], axis=3)
cfgs[0]["gpuid"] = -1
b2dcpu = np.sum(pmmc.mmclab(cfgs[0])["flux"], axis=3)

cfgs[1]["gpuid"] = 1
b2gpu = np.sum(pmmc.mmclab(cfgs[1])["flux"], axis=1)
cfgs[1]["gpuid"] = -1
b2cpu = np.sum(pmmc.mmclab(cfgs[1])["flux"], axis=1)

xi, yi = np.meshgrid(np.arange(0.5, 60), np.arange(0.5, 60))
plane = np.array([[0, 30.5, 0], [0, 30.5, 1], [1, 30.5, 0]])
interp = []
for b2 in (b2cpu, b2gpu):
    cutpos, cutvalue = i2m.qmeshcut(cfgs[1]["elem"][:, :4], cfgs[1]["node"], b2, plane)[:2]
    cutpos = np.asarray(cutpos)
    interp.append(griddata((cutpos[:, 0], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi)))
b2cpu_interpolated, b2gpu_interpolated = interp

## MCX-CL on the voxelized domain

cfg_mcx = {}
cfg_mcx["nphoton"] = cfgs[0]["nphoton"]
dim = 60
vx, vy, vz = np.meshgrid(*([np.arange(0.5, dim)] * 3), indexing="ij")
dist = (vx - 30) ** 2 + (vy - 30) ** 2 + (vz - 30) ** 2
vol = np.ones(vx.shape)
vol[dist < 625] = 2  # radius 25
vol[dist < 529] = 3  # radius 23
vol[dist < 100] = 4  # radius 10
cfg_mcx["vol"] = vol.astype("uint8")
cfg_mcx["srcpos"] = [30, 30.1, 0]
cfg_mcx["srcdir"] = [0, 0, 1]
cfg_mcx["issrcfrom0"] = 1
cfg_mcx["prop"] = [[0, 0, 1, 1], [0.02, 7.0, 0.89, 1.37], [0.0004, 0.009, 0.89, 1.37], [0.02, 9.0, 0.89, 1.37], [0.05, 0.0, 1.0, 1.37]]
cfg_mcx["tstart"] = 0
cfg_mcx["tend"] = 5e-9
cfg_mcx["tstep"] = 5e-10
cfg_mcx["autopilot"] = 1
cfg_mcx["gpuid"] = 1
cfg_mcx["isreflect"] = 1  # disable reflection at exterior boundaries
cfg_mcx["unitinmm"] = 1  # resolution
cfg_mcx["issavedet"] = 0  # mcxlabcl with one output
phimcx = pmcxcl.mcxlab(cfg_mcx)
b2mcx = np.sum(phimcx["flux"], axis=3)

## compare on y=30.5

tstep = cfgs[0]["tstep"]
clines = np.arange(-8, 0.01, 0.5)
plt.figure()
plt.contourf(xi, yi, np.log10(np.maximum(b2cpu_interpolated * tstep, 1e-8)), clines)
plt.contour(xi, yi, np.log10(np.abs(b2gpu_interpolated * tstep)), clines, colors="g", linestyles="--", linewidths=2)
plt.contour(xi, yi, np.log10(np.abs(b2dcpu[:-1, 30, :-1].T * tstep)), clines, colors="b", linestyles="--", linewidths=2)
plt.contour(xi, yi, np.log10(np.abs(b2dgpu[:-1, 30, :-1].T * tstep)), clines, colors="r", linewidths=1)
plt.contour(xi, yi, np.log10(np.abs(b2mcx[:, 30, :].T * cfg_mcx["tstep"])), clines, colors="m", linestyles="--")
t = np.linspace(0, 2 * np.pi, 200)
for r in (10, 23, 25):  # plotting at y=30.5, 0.5 mm from the center of the sphere
    rr = np.sqrt(r**2 - 0.5**2)
    plt.plot(rr * np.cos(t) + 30, rr * np.sin(t) + 30, "k--", linewidth=2)
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("B2 CPU (filled), B2 GPU (green), B2D CPU (blue), B2D GPU (red), MCX-CL (magenta)")
plt.show()
