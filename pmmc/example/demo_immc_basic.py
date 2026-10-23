"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_immc_basic.m

Implicit MMC (iMMC): thin vessels and membranes represented implicitly as an edge (edgeroi),
a node (noderoi) or a face (faceroi) region of interest, with radius/thickness 0.1, in a 1 mm
cube; iMMC runs on the CPU (gpuid=-1). This demo requires iso2mesh (pip install iso2mesh).

The same three cases as the benchmarks in the iMMC paper: https://doi.org/10.1364/BOE.411898
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## shared settings

cfg0 = {}
cfg0["nphoton"] = 1e6
cfg0["srcpos"] = [0.5, 0.5, 1]
cfg0["srcdir"] = [0, 0, -1]
cfg0["prop"] = [[0, 0, 1, 1], [0.0458, 35.6541, 0.9000, 1.3700], [23.0543, 9.3985, 0.9000, 1.3700]]
cfg0["tstart"] = 0
cfg0["tend"] = 5e-9
cfg0["tstep"] = 5e-9
cfg0["debuglevel"] = "TP"
cfg0["method"] = "grid"
cfg0["steps"] = [0.01, 0.01, 0.01]
cfg0["isreflect"] = 1
cfg0["gpuid"] = -1

## edge-based iMMC: a vessel along the edge between node 9 and 10

nbox, ebox = i2m.meshgrid6([0, 1], [0, 1], [0, 1])[:2]
fbox = i2m.volface(ebox)[0]
EPS = 0.001
nbox = np.vstack([nbox, [1 - EPS, 0.5, 0.5], [EPS, 0.5, 0.5]])  # insert new nodes (node 9 and 10)
fbox = [list(f) for f in fbox] + [[9, 9, 10]]  # insert new edge connected by node 9 and 10

cfg = copy.copy(cfg0)
cfg["node"], cfg["elem"] = i2m.s2m(nbox, fbox, 1, 100, "tetgen1.5", None, None, "-YY")[:2]
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
ne = cfg["elem"].shape[0]
edges = np.sort(np.asarray(i2m.meshedge(cfg["elem"][:, :4])), axis=1)
edgeid = np.reshape(np.all(edges == [9, 10], axis=1), (ne, 6), order="F")
cfg["edgeroi"] = np.zeros((ne, 6))
cfg["edgeroi"][edgeid] = 0.1
cfg["noderoi"] = np.zeros(cfg["node"].shape[0])
cfg["noderoi"][8] = 0.1
flux_eimmc = pmmc.mmclab(cfg)

## node-based iMMC: a sphere centered at node 9

nbox, ebox = i2m.meshgrid6([0, 1], [0, 1], [0, 1])[:2]
fbox = i2m.volface(ebox)[0]
nbox = np.vstack([nbox, [0.5, 0.5, 0.5]])  # insert new nodes (node 9)
cfg = copy.copy(cfg0)
cfg["node"], cfg["elem"] = i2m.s2m(nbox, [list(f) for f in fbox], 1, 100, "tetgen1.5", None, None, "-YY")[:2]
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["noderoi"] = np.zeros(cfg["node"].shape[0])
cfg["noderoi"][8] = 0.1
flux_nimmc = pmmc.mmclab(cfg)

## face-based iMMC: a membrane on the faces whose nodes are all beyond node 4

cfg = copy.copy(cfg0)
cfg["node"], cfg["elem"] = i2m.meshgrid6([0, 1], [0, 1], [0, 1])[:2]
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
facelist = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]  # nchoosek(1:4, 3)
allfaces = np.vstack([cfg["elem"][:, f] for f in facelist])
cfg["faceroi"] = np.reshape((np.sum(allfaces > 4, axis=1) == 3).astype(float), (cfg["elem"].shape[0], len(facelist)), order="F")
cfg["faceroi"][cfg["faceroi"] > 0] = 0.1
flux_fimmc = pmmc.mmclab(cfg)

## plot the fluence on the x=0.5 slice of the output grid

fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, (res, name) in zip(axs, ((flux_eimmc, "edge-based iMMC"), (flux_nimmc, "node-based iMMC"), (flux_fimmc, "face-based iMMC"))):
    im = ax.imshow(np.log10(np.rot90(np.abs(np.squeeze(res["flux"][49, :100, :100])))))
    ax.set_title(name)
    plt.colorbar(im, ax=ax)
plt.show()
