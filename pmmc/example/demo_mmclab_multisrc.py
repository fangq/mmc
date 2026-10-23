"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_mmclab_multisrc.m

Multiple sources in one simulation: cfg['srcpos'] and cfg['srcdir'] with one row per source
(the 4th column of srcpos is the source weight) and cfg['srcid']=-1 return the fluence of each
source separately. The first slot is then compared with a run of source 1 alone (srcid=1). This
demo requires iso2mesh (pip install iso2mesh).
"""

import copy
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshplot import plotcut

## prepare simulation input

node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 30], 6)

cfg = {}
cfg["nphoton"] = 1e7
cfg["node"] = node
cfg["elem"] = elem
cfg["elemprop"] = np.ones(elem.shape[0])
cfg["srcpos"] = [[20, 30, 0, 0.5], [40, 30, 0, 0.5]]
cfg["srcdir"] = [[0, 0, 1, 0], [0, 0, 1, 0]]
cfg["srcid"] = -1  # simulate both slots, separate outputs
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1, 0, 1.37]]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["method"] = "elem"  # mesh mode (Branchless Badouel)
cfg["basisorder"] = 1  # nodal fluence
cfg["debuglevel"] = "TP"

## run both sources together

print("\n-- Running multi-source forward (srcid=-1, 2 sources) --")
fluxall = pmmc.mmclab(cfg)["flux"]
print("  flux size: %s   nsrcslots=%d" % (str(fluxall.shape), fluxall.shape[-1]))

## run the first source alone with half of the photons

cfg1 = copy.copy(cfg)
cfg1["srcid"] = 1  # only the first slot (matches srcpos row 1)
cfg1["nphoton"] = cfg["nphoton"] / 2  # half the budget -> matches the per-source share above
print("\n-- Running srcid=1 only (slot 1, %.0e photons) --" % cfg1["nphoton"])
flux1 = pmmc.mmclab(cfg1)["flux"]

## plot and compare

slab = np.abs(fluxall[:, 0, :])  # [nn, 2]
fig, axs = plt.subplots(1, 3, figsize=(16, 4))
for k in range(2):
    plotcut(cfg["node"], cfg["elem"], np.log10(slab[:, k] + 1e-12), "y=30", axs[k], "source %d at [%g %g %g], weight %g" % ((k + 1,) + tuple(cfg["srcpos"][k])))
phi_srcid1 = np.abs(np.reshape(flux1, (flux1.shape[0], -1))[:, 0])
plotcut(cfg["node"], cfg["elem"], np.log10(phi_srcid1 + 1e-12), "y=30", axs[2], "srcid=1 alone")
fig.suptitle("Multi-source forward via cfg['srcpos'] = [...; ...]")

mask = slab[:, 0] > 1e-6
relerr = np.abs(phi_srcid1[mask] - slab[mask, 0]) / slab[mask, 0]
print("\nSlot-1 multi vs srcid=1 single-run agreement:")
print("  median rel.err: %.2f %%   95th pct: %.2f %%" % (100 * np.median(relerr), 100 * np.percentile(relerr, 95)))
plt.show()
