"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_optix_mmc.m

Compare the OptiX (hardware ray-tracing) and the OpenCL backends of MMC on a sphere in a box,
with dual-grid output (method='grid'). As in mmclab, compute='optix' runs on OpenCL when pmmc is
not built with OptiX. This demo requires iso2mesh (pip install iso2mesh).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## create the mesh: a sphere in a box

no_box, el_box = i2m.meshgrid6([0, 60], [0, 60], [0, 60])[:2]
fc_box = i2m.volface(el_box)[0]
no_sphere, fc_sphere = i2m.meshasphere([30, 30, 30], 25, 3, dovol=False)[:2]  # the surface, as [node, face] = meshasphere()
no, fc = i2m.mergemesh(no_box, fc_box, no_sphere, fc_sphere)[:2]
node, elem, face = i2m.surf2mesh(no, fc[:, :3], [0, 0, 0], [60, 60, 60], 1.0, 100, regions=np.array([[1, 1, 1], [30, 30, 30]]), cmdopt="-Y -A")

## simulation settings

cfg = {}
cfg["nphoton"] = 1e8
cfg["node"] = node
cfg["elem"] = elem
cfg["srcpos"] = [30, 30, 0.01]
cfg["srcdir"] = [0, 0, 1]
cfg["prop"] = [[0.000, 0, 1, 1], [0.005, 1, 0, 1.37], [0.010, 10, 0, 1]]  # box, sphere
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10
cfg["method"] = "grid"
cfg["outputtype"] = "fluence"
cfg["gpuid"] = 1

## run with OptiX and with OpenCL

cfg["compute"] = "optix"
fluenceoptix = np.sum(pmmc.mmclab(cfg)["flux"][:-1, :-1, :-1, :], axis=3)

cfg["compute"] = "opencl"
fluenceopencl = np.sum(pmmc.mmclab(cfg)["flux"][:-1, :-1, :-1, :], axis=3)

## compare on the x=30 slice

clines = np.arange(-25, -4.99, 2)
plt.figure()
plt.contourf(np.log(fluenceopencl[29, :, :].T), clines)
plt.contour(np.log(fluenceoptix[29, :, :].T), clines, colors="r", linestyles="--")
plt.gca().set_aspect("equal")
plt.title("log fluence at x=30, filled: OpenCL, red: OptiX")
plt.show()
