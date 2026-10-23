"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmclab_slit.m

A slit source (a line source emitting along a direction) outside of a 60x60x30 mm box; pmmc.mmclab
adds the source domain to the mesh. This demo requires iso2mesh (pip install iso2mesh).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshplot import plotcut

## prepare simulation input

cfg = {}
cfg["nphoton"] = 1e7
cfg["node"], face, cfg["elem"] = i2m.meshabox([0, 0, 0], [60, 60, 30], 6)
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [20, 30, -10]
cfg["srcdir"] = [0, 0.2, np.sqrt(1 - 0.04)]
cfg["prop"] = [[0, 0, 1, 1], [0.005, 0.1, 0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["debuglevel"] = "TP"
cfg["unitinmm"] = 1
cfg["method"] = "elem"
cfg["srctype"] = "slit"
cfg["srcparam1"] = [30, 0, 0, 0]
cfg["srcparam2"] = [0, 0, 0, 0]

## run the simulation

res = pmmc.mmclab(cfg)
flux = res["flux"]

## plot the fluence on two cut planes of the original mesh

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
logflux = np.log10(np.abs(flux[: cfg["node"].shape[0], 0]))
plotcut(cfg["node"], cfg["elem"], logflux, "x=40", axs[0], "log10 fluence at x=40")
plotcut(cfg["node"], cfg["elem"], logflux, "y=30", axs[1], "log10 fluence at y=30")
plt.show()
