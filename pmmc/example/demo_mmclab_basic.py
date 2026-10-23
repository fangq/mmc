"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmclab_basic.m

A basic MMC simulation: a pencil beam into a 60x60x30 mm box, saving the volumetric fluence and
the diffuse reflectance on the surface. This demo requires iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## prepare simulation input

cfg = {}
cfg["nphoton"] = 1e6
cfg["node"], face, cfg["elem"] = i2m.meshabox([0, 0, 0], [60, 60, 30], 6)
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [30, 30, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1, 0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["debuglevel"] = "TP"
cfg["issaveref"] = 1  # in addition to volumetric fluence, also save surface diffuse reflectance
cfg["method"] = "elem"

## run the simulation

flux = pmmc.mmclab(cfg)

## plotting the result

# fluence on the y=30 plane
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
cutpos, cutval = i2m.qmeshcut(cfg["elem"][:, :4], cfg["node"], np.log10(np.abs(flux["flux"][: cfg["node"].shape[0], 0])), "y=30")[:2]
cutpos = np.asarray(cutpos)
tc = axs[0].tricontourf(cutpos[:, 0], cutpos[:, 2], np.ravel(cutval), 40)
axs[0].invert_yaxis()
axs[0].set_title("log10 fluence at y=30")
fig.colorbar(tc, ax=axs[0])

# diffuse reflectance on the surface faces, shown at the face centroids of the bottom (z=0)
if cfg.get("issaveref") == 1:
    faces = i2m.faceneighbors(cfg["elem"][:, :4], "rowmajor")
    faces = faces[0] if isinstance(faces, tuple) else faces
    centroid = cfg["node"][faces - 1].mean(axis=1)
    bottom = centroid[:, 2] < 1e-6
    sc = axs[1].scatter(centroid[bottom, 0], centroid[bottom, 1], c=np.log10(flux["dref"][bottom, 0]), s=8)
    axs[1].set_title("log10 diffuse reflectance on z=0")
    axs[1].set_aspect("equal")
    fig.colorbar(sc, ax=axs[1])

plt.show()
