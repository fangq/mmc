"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmclab_mesh_adjoint.m

Mesh-mode (nodal) adjoint Jacobians: with outputtype='adjoint_mua_d', one simulation launches
photons from the source and from each detector, and returns the mua Jacobian (res['jmua']) and
the diffusion-coefficient Jacobian (res['jd']) of every source-detector pair, [nn, Ns*Nd]. This
demo requires iso2mesh (pip install iso2mesh).
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
cfg["detpos"] = [[20, 30, 0, 1.5], [40, 30, 0, 1.5], [45, 30, 0, 1.5]]
cfg["detdir"] = [[0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0]]
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1, 0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["method"] = "elem"  # mesh mode (Branchless Badouel)
cfg["basisorder"] = 1  # nodal fluence required for J_D
cfg["outputtype"] = "adjoint_mua_d"
cfg["adjointmode"] = 0  # 0 = full FEM (default), 1 = nodal approximation
cfg["debuglevel"] = "TP"

## run the simulation

print("\n-- Running CW mesh-mode adjoint (Ns=1, Nd=3 -> 3 J_mua + 3 J_D) --")
flux = pmmc.mmclab(cfg)
print("  flux size: %s" % str(flux["flux"].shape))
print("  jmua size: %s   range: [%.2e, %.2e]" % (str(flux["jmua"].shape), flux["jmua"].min(), flux["jmua"].max()))
print("  jd   size: %s   range: [%.2e, %.2e]" % (str(flux["jd"].shape), flux["jd"].min(), flux["jd"].max()))

## plot the Jacobians of each source-detector pair

npair = flux["jmua"].shape[1]
fig, axs = plt.subplots(2, npair, figsize=(15, 7))
for k in range(npair):
    sd = cfg["detpos"][k][0] - cfg["srcpos"][0]
    plotcut(cfg["node"], cfg["elem"], np.log10(np.abs(flux["jmua"][:, k]) + 1e-12), "y=30", axs[0, k], "log10|J_mua|, S-D %g mm" % sd)
    plotcut(cfg["node"], cfg["elem"], np.log10(np.abs(flux["jd"][:, k]) + 1e-12), "y=30", axs[1, k], "log10|J_D|, S-D %g mm" % sd)
fig.suptitle("Mesh-mode adjoint (full FEM): J_mua (top) and J_D (bottom) for 3 S-D pairs")

## consistency: J_mua(S, D1) vs phi_S * phi_D1

phi_src = np.reshape(flux["flux"], (flux["flux"].shape[0], -1))[:, 0]  # slot 1 = forward S
phi_d1 = np.reshape(flux["flux"], (flux["flux"].shape[0], -1))[:, 1]  # slot 2 = first detector adjoint
prod_sd = phi_src * phi_d1
mask = (np.abs(prod_sd) > 1e-12) & (np.abs(flux["jmua"][:, 0]) > 1e-12)
cc = np.corrcoef(np.log10(np.abs(prod_sd[mask])), np.log10(np.abs(flux["jmua"][mask, 0])))
print("\nConsistency: corrcoef( J_mua(S,D_1) , phi_S * phi_D_1 ):")
print("  log10-log10 corr coef = %.4f  (expected close to 1.0)" % cc[0, 1])
plt.show()
