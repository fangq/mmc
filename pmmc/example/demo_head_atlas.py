"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_head_atlas.m

Dual-grid MMC (DMMC, method='grid') in a brain atlas mesh (head_atlas.mat in
mmclab/example), with a pencil beam at the EEG 10-5 landmark C4h. This demo requires iso2mesh
(pip install iso2mesh) and scipy.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.io import loadmat
import pmmc

## load the mesh

atlas = loadmat(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../mmclab/example/head_atlas.mat"))
node, elem, prop = np.double(atlas["node"]), atlas["elem"].astype(int), np.ravel(atlas["prop"]).astype(int)

cfg = {}
cfg["nphoton"] = 1e7
cfg["node"] = node
cfg["elem"] = elem
cfg["elemprop"] = prop
cfg["prop"] = [[0, 0, 1, 1], [0, 0, 1, 1], [0.019, 7.8, 0.89, 1.37], [0.019, 7.8, 0.89, 1.37], [0.0004, 0.009, 0.89, 1.37], [0.02, 9.0, 0.89, 1.37], [0.08, 40.9, 0.84, 1.37]]
cfg["srctype"] = "pencil"
cfg["srcdir"] = np.array([-0.5086, -0.1822, -0.8415])  # inward-pointing source
cfg["srcpos"] = np.array([133.5370, 90.1988, 200.0700])  # pencil beam source placed at EEG 10-5 landmark:"C4h"
cfg["srcpos"] = cfg["srcpos"] + 0.001 * cfg["srcdir"]  # ensure source is inside the mesh
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["isreflect"] = 1
cfg["debuglevel"] = "TP"
cfg["method"] = "grid"  # DMMC mode

## run the simulation

flux = pmmc.mmclab(cfg)
CWfluence = np.sum(flux["flux"] * cfg["tstep"], axis=3)

## plot the CW fluence on the y=90.5 plane, with the tissue boundaries

y_plane = 90.5
DMMC_AABB_origin = node.min(axis=0)  # spatial offset of the DMMC AABB
DMMC_voxel_origin = DMMC_AABB_origin + 0.5  # centroid of the 1st voxel for DMMC AABB
dims = CWfluence.shape
CWfluence_DMMC = CWfluence[:, int(np.ceil(y_plane - DMMC_AABB_origin[1])) - 1, :]
xx, zz = np.meshgrid(DMMC_voxel_origin[0] + np.arange(dims[0]), DMMC_voxel_origin[2] + np.arange(dims[2]))

plt.figure(figsize=(8, 7))
plt.contourf(xx, zz, np.log10(np.abs(CWfluence_DMMC.T)), np.arange(-20, 0.01, 0.5))
plt.colorbar()

plane = np.array([[0, y_plane, 0], [0, y_plane, 1], [1, y_plane, 0]])
label = np.unique(prop)
face = np.vstack([i2m.volface(elem[prop == label[i], :4])[0] for i in range(0, len(label), 2)])
cutpos, cutvalue, cutedges = i2m.qmeshcut(face[:, :3], node, node[:, 0], plane)[:3]
cutpos = np.asarray(cutpos)
plt.plot(cutpos[:, 0], cutpos[:, 2], "k.", markersize=1)
plt.plot(cfg["srcpos"][0], cfg["srcpos"][2], "o", markeredgecolor="r", markerfacecolor="r", markersize=10, label="source")
plt.gca().set_aspect("equal")
plt.title("DMMC CW fluence (log10) at y=%g" % y_plane)
plt.legend(loc="upper right")
plt.show()
