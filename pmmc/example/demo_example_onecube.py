"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_example_onecube.m

A 10x10x10 mm cube split into 6 tetrahedra, simulated with 100 photons, once with the photon
movement debug flag ('M') and once with the accumulation flag ('A'). This demo requires iso2mesh
(pip install iso2mesh).
"""

import numpy as np
import iso2mesh as i2m
import pmmc

## prepare simulation input

cfg = {}
cfg["node"], cfg["elem"] = i2m.meshgrid6([0, 10], [0, 10], [0, 10])[:2]  # genT6mesh in MATLAB
cfg["elem"] = cfg["elem"][np.lexsort(cfg["elem"].T[::-1])]  # sortrows
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1.0101010101, 0.01, 1.0]]
cfg["srcpos"] = [2, 8, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["seed"] = 1648335518
cfg["nphoton"] = 100
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-10

## run the simulation with photon movement and accumulation debug output

cfg["debuglevel"] = "M"
flux = pmmc.mmclab(cfg)

cfg["debuglevel"] = "A"
flux = pmmc.mmclab(cfg)
