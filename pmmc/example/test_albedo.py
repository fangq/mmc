"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/test_albedo.m

Compare the two photon weight schemes of MMC: albedo-weight (mcmethod=1, as MCML) and
micro-Beer-Lambert (mcmethod=0, as MCX), on the fluence and on the detected photons. This demo
requires iso2mesh (pip install iso2mesh).
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.interpolate import griddata
import pmmc

## prepare simulation input

cfg = {}
cfg["nphoton"] = 1e7
cfg["seed"] = -1
cfg["node"], cfg["elem"] = i2m.meshgrid6(np.arange(0, 61, 2), np.arange(0, 61, 2), np.arange(0, 61, 2))[:2]  # genT6mesh in MATLAB
cfg["elemprop"] = np.ones(cfg["elem"].shape[0])
cfg["srcpos"] = [30.1, 30.2, 0]
cfg["srcdir"] = [0, 0, 1]
cfg["tstart"] = 0
cfg["tend"] = 5e-9
cfg["tstep"] = 5e-9
cfg["prop"] = [[0, 0, 1, 1], [0.1, 1.0, 0.0, 1.0]]
cfg["debuglevel"] = "TP"
cfg["isreflect"] = 0
cfg["detpos"] = [[30.1, 25.2, 0.0, 1.0]]
cfg["issaveexit"] = 1  # save detected photon exit position and angles
cfg["issaveseed"] = 0  # mmclab with 2 outputs does not save the seeds
cfg["issavedet"] = 1  # mmclab saves the detected photons when its 2nd output is requested
cfg["isnormalized"] = 1
cfg["method"] = "elem"

## run with albedo-weight (mcmethod=1) and micro-Beer-Lambert (mcmethod=0)

cfg["mcmethod"] = 1
newcfg = pmmc.mmclab(cfg, "prep")  # preprocessing of the mesh to get the missing fields
res1 = pmmc.mmclab(newcfg)
cube1, detp1 = res1["flux"], res1["detp"]

cfg["mcmethod"] = 0
newcfg = pmmc.mmclab(cfg, "prep")
res2 = pmmc.mmclab(newcfg)
cube2, detp2 = res2["flux"], res2["detp"]

## compare the fluence on x=30

xi, yi = np.meshgrid(np.arange(61), np.arange(61))
plane = np.array([[30, 0, 0], [30, 0, 1], [30, 1, 0]])
vis = []
for cube in (cube1, cube2):
    cutpos, cutvalue = i2m.qmeshcut(cfg["elem"][:, :4], cfg["node"], np.ravel(cube[:, 0]), plane)[:2]
    cutpos = np.asarray(cutpos)
    vis.append(griddata((cutpos[:, 1], cutpos[:, 2]), np.ravel(cutvalue), (xi, yi)))

plt.figure()
plt.contour(np.log10(np.abs(vis[0])), np.arange(1, 8.01, 0.5), colors="r", linestyles=":", linewidths=1.5)
plt.contour(np.log10(np.abs(vis[1])), np.arange(1, 8.01, 0.5), colors="b", linewidths=1.5)
plt.gca().invert_yaxis()
plt.gca().set_aspect("equal")
plt.title("fluence at x=30, red: albedo-weight, blue: micro-Beer-Lambert")

## compare the detected photon weights, mean scattering counts and partial paths

mua, mus = 0.1, 1
mut = mua + mus
w1 = detp1["w0"] * (mus / mut) ** detp1["nscat"].astype(float).ravel()
w2 = detp2["w0"] * np.exp(-mua * detp2["ppath"].astype(float).ravel())
print("relative difference of the total detected weight: %g" % (abs(np.sum(w1) - np.sum(w2)) / np.sum(w2)))

c0, n = 3e11, 1
cc = c0 / n
nTG, gate_width = 50, 2e-11
tpsf1, tpsf2 = np.zeros(nTG), np.zeros(nTG)
for ppath, w, tpsf in ((detp1["ppath"].ravel(), w1, tpsf1), (detp2["ppath"].ravel(), w2, tpsf2)):
    tg = np.floor(ppath / cc / gate_width).astype(int)
    keep = tg < nTG
    np.add.at(tpsf, tg[keep], w[keep])

plt.figure()
plt.plot(tpsf1, "r", tpsf2, "b")
plt.legend(["albedo-weight", "micro-Beer-Lambert"])
plt.title("TPSF of the detected photons")

p1_ave = np.sum(detp1["nscat"].astype(float).ravel() * w1) / len(w1)
p2_ave = np.sum(detp2["nscat"].astype(float).ravel() * w2) / len(w2)
L1_ave = np.sum(detp1["ppath"].astype(float).ravel() * w1) / len(w1)
L2_ave = np.sum(detp2["ppath"].astype(float).ravel() * w2) / len(w2)
print("mean scattering count: %g vs %g, mean path length: %g vs %g" % (p1_ave, p2_ave, L1_ave, L2_ave))
plt.show()
