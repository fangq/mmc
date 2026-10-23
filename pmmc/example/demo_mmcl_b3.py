"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmcl_b3.m

Benchmark B3 of the MMC paper: an adult brain atlas mesh (Collins atlas,
MMC_Collins_Atlas_Mesh_Version_2L.mat in mmclab/example), simulated on the GPU (gpuid=1) and on
the CPU (gpuid=-1). This demo requires iso2mesh (pip install iso2mesh) and scipy.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.io import loadmat
from scipy.interpolate import griddata
import pmmc

## load the mesh

matfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../mmclab/example/MMC_Collins_Atlas_Mesh_Version_2L.mat")
if not os.path.isfile(matfile):
    raise FileNotFoundError(
        "please download MMC_Collins_Atlas_Mesh_Version_2L.tar.gz from https://sourceforge.net/projects/mcx/files/mmc/ "
        "and extract the .mat file to mmclab/example"
    )
atlas = loadmat(matfile)
node, elem, face = atlas["node"].astype(float), atlas["elem"].astype(int), atlas["face"].astype(int)

cfg = {}
cfg["node"] = node
cfg["elem"] = elem[:, :4]
cfg["elemprop"] = elem[:, 4]
cfg["tstart"] = 0
cfg["tend"] = 5e-09
cfg["tstep"] = 2e-10
cfg["srcpos"] = [75, 67.38, 167.5]
srcdir = np.array([0.1636, 0.4569, -0.8743])
cfg["srcdir"] = srcdir / np.linalg.norm(srcdir)
cfg["detpos"] = [[75.0000, 77.1900, 170.3000, 3.0000], [75.0000, 89.0000, 171.6000, 3.0000], [75.0000, 97.6700, 172.4000, 3.0000], [75.0000, 102.4000, 172.0000, 3.0000]]
cfg["prop"] = [
    [0, 0, 1.0000, 1.0000],  # background/air
    [0.0190, 7.8182, 0.8900, 1.3700],  # scalp & skull are both labeled using 1
    [0.0004, 0.0090, 0.8900, 1.3700],  # csf
    [0.0200, 9.0000, 0.8900, 1.3700],  # gray matters
    [0.0800, 40.9000, 0.8400, 1.3700],  # white matters
    [0, 0, 1.0000, 1.0000],  # air pockets
]
cfg["seed"] = 29012392
cfg["nphoton"] = 1e8
cfg["method"] = "elem"

## run on the GPU and on the CPU

cfg["gpuid"] = 1
b3gpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=1) * cfg["tstep"]
cfg["gpuid"] = -1
b3cpu = np.sum(pmmc.mmclab(cfg)["flux"], axis=1) * cfg["tstep"]

## compare on the x=74.5 plane

lo, hi = node.min(axis=0), node.max(axis=0)
bb_size = np.floor(hi - lo).astype(int)
yi, zi = np.meshgrid(lo[1] + 0.5 + np.arange(bb_size[1]), lo[2] + 0.5 + np.arange(bb_size[2]), indexing="ij")
plane = np.array([[74.5, 0, 0], [74.5, 0, 1], [74.5, 1, 0]])
interp = []
for b3 in (b3cpu, b3gpu):
    cutpos, cutvalue = i2m.qmeshcut(elem[:, :4], node, b3, plane)[:2]
    cutpos = np.asarray(cutpos)
    v = griddata((cutpos[:, 1], cutpos[:, 2]), np.ravel(cutvalue), (yi, zi))
    v[v == 0] = np.nan
    interp.append(v)

clines = np.arange(-10, 0.01, 0.5)
plt.figure(figsize=(8, 7))
plt.contourf(yi, zi, np.log10(interp[0]), clines)
plt.contour(yi, zi, np.log10(interp[1]), clines, colors="g", linestyles="--", linewidths=2)
cutpos = np.asarray(i2m.qmeshcut(face[:, :3], node, node[:, 0], np.array([[74.5, 0, 0], [74.5, 0, 1], [74.5, 1, 0]]))[0])
plt.plot(cutpos[:, 1], cutpos[:, 2], "k.", markersize=1)
plt.gca().set_aspect("equal")
plt.title("B3 CW fluence (log10), filled: CPU, green: GPU")
plt.show()
