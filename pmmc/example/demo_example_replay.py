"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_example_replay.m

Replay the photons detected by one detector to compute its Jacobian: the initial run saves the
detected photons and their RNG seeds, and the replay runs with outputtype 'wl' and 'wp' rerun
only those photons. This demo requires iso2mesh (pip install iso2mesh).
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
import pmmc

## run the initial simulation

cfg = {}
cfg["nphoton"] = 1e6
cfg["seed"] = 1648335518
cfg["node"], cfg["elem"] = i2m.meshgrid6(np.arange(0, 61, 2), np.arange(0, 61, 2), np.arange(0, 61, 2))[:2]  # genT6mesh in MATLAB
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [30.1, 30.2, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["prop"] = [[0, 0, 1, 1], [0.005, 1.0, 0.01, 1.0]]
cfg["debuglevel"] = "TP"
cfg["method"] = "elem"
cfg["isreflect"] = 0
cfg["detpos"] = [[30.0, 20.0, 0.0, 1.0], [30.0, 40.0, 0.0, 1.0], [20.0, 30.0, 0.0, 1.0], [40.0, 30.0, 0.0, 1.0]]
cfg["issaveexit"] = 1  # save detected photon exit position and angles
cfg["issaveseed"] = 1  # save detected photon seeds to replay later
cfg["issavedet"] = 1  # mmclab saves the detected photons when its 2nd output is requested

newcfg = pmmc.mmclab(cfg, "prep")  # preprocessing of the mesh to get the missing fields
res = pmmc.mmclab(newcfg)  # initial simulation
cube, detp, seeds = res["flux"], res["detp"], res["seeds"]

## replay the photons detected by det#1, mua Jacobian ('wl')

newcfg = copy.copy(newcfg)
newcfg["replaydet"] = 1  # replay photons detected by det#1
newcfg["seed"] = seeds
newcfg["detphotons"] = detp["data"]
newcfg["isnormalized"] = 0
newcfg["outputtype"] = "wl"  # replay and get wl
newcfg["issaveseed"] = 0  # mmclab with 2 outputs does not save the seeds
res2 = pmmc.mmclab(newcfg)
cube2, detp2 = res2["flux"], res2["detp"]


def isreplayed(det0, det1):
    """every detected photon of the baseline run is found again in the replay"""
    a = {tuple(r) for r in np.round(det1["data"].T * 1e10) * 1e-10}
    return all(tuple(r) in a for r in np.round(det0.T * 1e10) * 1e-10)


print("replay is successful :-)" if isreplayed(detp["data"], detp2) else "replay failed :-(")


def plotcube(data, title):
    fig, axs = plt.subplots(1, 2, figsize=(11, 5))
    v = np.log10(np.abs(np.ravel(data[:, 0])) + 1e-30)
    for ax, plane, (i, j) in ((axs[0], "z=0.2", (0, 1)), (axs[1], "x=30", (1, 2))):
        cutpos, cutval = i2m.qmeshcut(newcfg["elem"][:, :4], newcfg["node"], v, plane)[:2]
        cutpos = np.asarray(cutpos)
        tc = ax.tricontourf(cutpos[:, i], cutpos[:, j], np.ravel(cutval), 40)
        ax.set_title(title + " at " + plane)
        ax.set_aspect("equal")
        fig.colorbar(tc, ax=ax)


plotcube(cube2, "log10 replay wl")

## replay det#1 again for the mus Jacobian ('wp')

newcfg["outputtype"] = "wp"  # replay and get wp
res3 = pmmc.mmclab(newcfg)
cube3, detp3 = res3["flux"], res3["detp"]

det1 = detp["data"][:, detp["data"][0, :] == newcfg["replaydet"]]
print("replay is successful :-)" if isreplayed(det1, detp3) else "replay failed :-(")

plotcube(cube3, "log10 replay wp")
plt.show()
