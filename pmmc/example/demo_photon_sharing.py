"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_photon_sharing.m

Photon sharing: one simulation computes the fluence of 5 illumination patterns at once. pmmc.mmclab
adds the source domain of the pattern source to the mesh. This demo requires iso2mesh
(pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## prepare simulation input

cfg = {}
cfg["nphoton"] = 3e6
cfg["node"], face, cfg["elem"] = i2m.meshabox([0, 0, 0], [60, 60, 20], 2, 2)
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [10, 10, -2]
cfg["srcdir"] = [0, 0, 1]
cfg["srctype"] = "pattern"
cfg["srcparam1"] = [40.0, 0.0, 0.0, 40]
cfg["srcparam2"] = [0.0, 40.0, 0.0, 40]
cfg["prop"] = [[0, 0, 1, 1], [0.01, 10, 0.9, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 1e-10
cfg["gpuid"] = 1
cfg["method"] = "elem"
cfg["debuglevel"] = "TP"
cfg["issavedet"] = 1  # mmclab saves the detected photons when its 2nd output is requested

pat = np.zeros([40, 40, 5])
pat[:, :, 0] = 1  # full field illumination pattern
pat[0:20, :, 1] = 1  # pattern with left half bright
pat[:, 0:20, 2] = 1  # pattern with top half bright
pat[10:30, 10:30, 3] = 1  # pattern with bright square in the middle
pat[15:25, :, 4] = 1
pat[:, 15:25, 4] = 1  # pattern with a bright cross

cfg["srcpattern"] = np.transpose(pat, [2, 0, 1])  # [number of patterns, Nx, Ny]

## run the simulation

res = pmmc.mmclab(cfg)
flux = res["flux"]  # [number of patterns, number of nodes of the extended mesh, time gates]

## plot the results

node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 20], 2, 2)
data = flux[:, : node.shape[0], :]  # the nodes of the original mesh come first
cwdata = np.sum(data, axis=2)

fig, axs = plt.subplots(2, 3, figsize=(13, 8))
bottom = node[:, 2] < 1e-6
for k, ax in zip(range(5), [axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 1], axs[1, 2]]):
    sc = ax.tricontourf(node[bottom, 0], node[bottom, 1], cwdata[k, bottom], 40)
    ax.set_title("pattern %d" % (k + 1))
    ax.set_aspect("equal")
    ax.axis("off")
axs[1, 0].axis("off")
plt.show()
