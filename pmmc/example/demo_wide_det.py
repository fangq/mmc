"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_wide_det.m

A wide-field (40x40 mm planar) detector below a two-layer slab: the detected photons are binned
into a 40x40 image per time gate (issaveexit=2). This demo requires iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## create the two-layer mesh and add the detector domain

node, face, c0 = i2m.latticegrid([0, 60], [0, 60], [0, 5, 10])
c0 = np.hstack([np.asarray(c0, dtype=float)[:, :3], np.array([[2], [3]])])  # maximum element size of the bottom (label 1) and top (label 2) layers
node, elem = i2m.surf2mesh(node, face, None, None, 1, None, regions=c0)[:2]

detdef = {"srctype": "planar", "srcpos": [10, 10, -1], "srcdir": [0, 0, 1], "srcparam1": [40, 0, 0, 40], "srcparam2": [0, 40, 0, 40]}

cfg = {}
cfg["node"], cfg["elem"] = pmmc.mmcadddet(node, elem, detdef)
cfg["elemprop"] = cfg["elem"][:, 4]
cfg["elem"] = cfg["elem"][:, :4]

## define other properties

cfg["prop"] = [[0, 0, 1, 1], [0.01, 1.0, 0.01, 1.37], [0.05, 10.0, 0.9, 1.37]]
cfg["srcpos"] = [25.0, 35.0, 10.0]
cfg["e0"] = i2m.tsearchn(cfg["node"], cfg["elem"][:, :4], np.array([cfg["srcpos"]], dtype=float))[0][0]
cfg["srcdir"] = [0, 0, -1]
cfg["srctype"] = "pencil"
cfg["srcparam1"] = [0, 0, 0, 0]
cfg["srcparam2"] = [0, 0, 0, 0]
cfg["detpos"] = [[10.0, 10.0, -1.0, 0]]
cfg["detparam1"] = [40.0, 0.0, 0.0, 40]
cfg["detparam2"] = [0.0, 40.0, 0.0, 40]
cfg["tstart"] = 0
cfg["tend"] = 2e-9
cfg["tstep"] = 4e-11
cfg["nphoton"] = 1e6
cfg["seed"] = 12345678
cfg["debuglevel"] = "TP"
cfg["issaveexit"] = 2
cfg["issavedet"] = 1  # mmclab saves the detected photons when its 2nd output is requested

## run the simulation

res = pmmc.mmclab(cfg)
detp = res["detp"]  # with issaveexit=2, a [40, 40, time gates] image of the detector

## plot the result

plt.figure()
plt.imshow(np.sum(detp, axis=2).T, origin="upper")
plt.title("total detected weight on the wide-field detector")
plt.axis("equal")
plt.colorbar()
plt.show()
