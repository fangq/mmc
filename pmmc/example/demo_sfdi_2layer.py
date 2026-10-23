"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_sfdi_2layer.m

Spatial-frequency domain imaging (SFDI): a Fourier pattern source above a layered slab (3 layers
by default, set layercount = 2 for 2 layers). The source domain is added to the mesh explicitly
with mmcsrcdomain/mmcaddsrc. This demo requires iso2mesh (pip install iso2mesh).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshplot import plotcut

## create the layered mesh

layercount = 3

cfg = {}
if layercount == 2:
    node, face, c0 = i2m.latticegrid([0, 60], [0, 60], [0, 25, 30])
    sizes = [[2], [1]]  # maximum element size for bottom (label 1) and top (label 2) layers
    cfg["prop"] = [[0, 0, 1, 1], [0.02, 9.0, 0.89, 1.37], [0.019, 7.8, 0.89, 1.37]]
else:
    node, face, c0 = i2m.latticegrid([0, 60], [0, 60], [0, 20, 25, 30])  # if you like a 3-layer model
    sizes = [[2], [2], [1]]
    cfg["prop"] = [[0, 0, 1, 1], [0.02, 9.0, 0.89, 1.37], [0.0004, 0.009, 0.89, 1.37], [0.019, 7.8, 0.89, 1.37]]
c0 = np.hstack([np.asarray(c0, dtype=float)[:, :3], np.array(sizes)])

cfg["node"], cfg["elem"] = i2m.surf2mesh(node, face, None, None, 1, None, regions=c0)[:2]
cfg["elemprop"] = cfg["elem"][:, 4]
cfg["elem"] = cfg["elem"][:, :4]

## define the SFDI source and add its domain to the mesh

cfg["srctype"] = "fourier"  # define an SFDI source
cfg["srcpos"] = [10, 10, 35]  # one corner of the illumination area
kx = 3  # wave number in the x-dir
ky = 0  # wave number in the y-dir
xphase = np.pi / 3  # phase offset in the x-dir, must < 2pi
yphase = 0  # phase offset in the y-dir, must < 2pi
cfg["srcparam1"] = [40, 0, 0, kx + xphase / (2 * np.pi)]  # kx is k-number in x direction
cfg["srcparam2"] = [0, 40, 0, ky + yphase / (2 * np.pi)]
cfg["srcdir"] = [0, 0, -1]

srcdef = {k: cfg[k] for k in ("srctype", "srcpos", "srcdir", "srcparam1", "srcparam2")}
sdom = pmmc.mmcsrcdomain(srcdef, [cfg["node"].min(axis=0), cfg["node"].max(axis=0)])[0]
cfg["node"], cfg["elem"] = pmmc.mmcaddsrc(cfg["node"], np.hstack([cfg["elem"], cfg["elemprop"][:, None]]), sdom)
cfg["elemprop"] = cfg["elem"][:, 4]
cfg["elem"] = cfg["elem"][:, :4]

## other simulation settings

cfg["nphoton"] = 3e6
cfg["seed"] = 1648335518
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 1e-10
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 1
cfg["detpos"] = [[30, 30, 0, 2]]  # detector position
cfg["method"] = "elem"

## run the simulation

layer = pmmc.mmclab(cfg)["flux"]
layercw = np.sum(layer, axis=1)

## plot the CW fluence in the tissue

tissue = cfg["elem"][cfg["elemprop"] > 0]
planes = ["y=30", "z=27"] + (["z=22"] if layercount != 2 else [])
fig, axs = plt.subplots(1, len(planes), figsize=(5 * len(planes), 4))
for ax, plane in zip(axs, planes):
    plotcut(cfg["node"], tissue, np.log10(layercw), plane, ax, "log10 CW fluence at " + plane)
plt.show()
