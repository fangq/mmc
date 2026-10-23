"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmclab_mesh_rf_adjoint.m

Frequency-domain (100 MHz) mesh-mode adjoint mua Jacobian: with cfg['omega'] > 0 and
outputtype='adjoint', res['jmua'] and res['flux'] are complex. This demo requires iso2mesh
(pip install iso2mesh).
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
cfg["node"], face, cfg["elem"] = i2m.meshabox([0, 0, 0], [60, 60, 30], 4)
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["nphoton"] = 3e7
cfg["srcpos"] = [30, 30, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["detpos"] = [[20, 30, 0, 1.5], [40, 30, 0, 1.5]]
cfg["detdir"] = [[0, 0, 1, 0], [0, 0, 1, 0]]
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1, 0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["method"] = "elem"
cfg["basisorder"] = 1  # nodal fluence required for adjoint
cfg["omega"] = 2 * np.pi * 100e6
cfg["outputtype"] = "adjoint"  # res['jmua'] only (complex)
cfg["adjointmode"] = 0  # 0 = full FEM
cfg["debuglevel"] = "TP"

## run the simulation

print("\n-- Running RF mesh-mode adjoint at %.0f MHz --" % (cfg["omega"] / (2 * np.pi * 1e6)))
flux = pmmc.mmclab(cfg)
jmua = flux["jmua"]
print("  jmua complex: %d   size: %s" % (np.iscomplexobj(jmua), str(jmua.shape)))
print("  |jmua| range: [%.2e, %.2e]" % (np.abs(jmua).min(), np.abs(jmua).max()))
print("  flux complex: %d   size: %s" % (np.iscomplexobj(flux["flux"]), str(flux["flux"].shape)))

## plot amplitude and phase of each pair

npair = jmua.shape[1]
fig, axs = plt.subplots(2, npair, figsize=(12, 7))
for k in range(npair):
    plotcut(cfg["node"], cfg["elem"], np.log10(np.abs(jmua[:, k]) + 1e-12), "y=30", axs[0, k], "log10|J_mua| pair %d" % (k + 1))
    plotcut(cfg["node"], cfg["elem"], np.angle(jmua[:, k]) * 180 / np.pi, "y=30", axs[1, k], "angle(J_mua) (deg) pair %d" % (k + 1))
fig.suptitle("Mesh-mode RF (100 MHz) adjoint Jacobian J_mua: amplitude and phase")

## consistency with phi_S * phi_D1

f2 = np.reshape(flux["flux"], (flux["flux"].shape[0], -1))
prod_sd = f2[:, 0] * f2[:, 1]
mask = (np.abs(prod_sd) > 1e-12) & (np.abs(jmua[:, 0]) > 1e-12)
if np.any(mask):
    cc_amp = np.corrcoef(np.log10(np.abs(prod_sd[mask])), np.log10(np.abs(jmua[mask, 0])))
    print("\nConsistency (log|.| corrcoef):  |J_mua(S,D_1)| vs |phi_S * phi_D_1| = %.4f" % cc_amp[0, 1])
    dphase = np.mod(np.angle(jmua[mask, 0]) - np.angle(prod_sd[mask]) + np.pi, 2 * np.pi) - np.pi
    print("Phase difference (rad):  median %.4f   max |.| %.4f" % (np.median(dphase), np.max(np.abs(dphase))))
plt.show()
