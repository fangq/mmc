"""
pmmc - Mesh-based Monte Carlo (MMC) for Python, ported from mmclab/example/demo_immc_vessel.m

Edge-based implicit MMC (iMMC) in a vessel network (vessel.mat in mmclab/example): every vessel
segment is a mesh edge with a radius (edgeroi), so thin vessels need no dense mesh. iMMC runs on
the CPU (gpuid=-1). This demo requires iso2mesh (pip install iso2mesh) and scipy.
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m
from scipy.io import loadmat
import pmmc

EDGEORDER = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]  # edge order in an element
FACEORDER = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]  # face order in an element


def edge2radius(edge, noder):
    """the radius of each vessel edge (the larger of its two nodes), keyed by the sorted node pair"""
    edge = np.sort(edge, axis=1)
    return {tuple(e): max(noder[e[0] - 1], noder[e[1] - 1]) for e in edge}


def splitelem(node, elem, eroi, e2r, index):
    """split element index into 4 at its centroid, so that it has at most 2 vessel edges"""
    selem = elem[index, :4]
    cindex = node.shape[0] + 1
    node = np.vstack([node, node[selem - 1].mean(axis=0)])  # put the centroid in the end
    newelem, newroi = [], []
    for f in FACEORDER:
        e = [selem[f[0]], selem[f[1]], selem[f[2]], cindex]
        newelem.append(e)
        newroi.append([e2r.get(tuple(sorted((e[a], e[b]))), 0.0) for a, b in EDGEORDER])
    elem = np.vstack([elem[:index], newelem, elem[index + 1 :]])
    eroi = np.vstack([eroi[:index], newroi, eroi[index + 1 :]])
    return node, elem, eroi


def edgeroi(node, elem, e2r):
    """the edge-ROI radii of each element, splitting elements with more than 2 vessel edges"""
    eroi = np.zeros((elem.shape[0], 6))
    i = 0
    while i < elem.shape[0]:
        enum = 0  # number of edges being labeled
        for j, (a, b) in enumerate(EDGEORDER):
            key = tuple(sorted((elem[i, a], elem[i, b])))
            if key in e2r:
                eroi[i, j] = e2r[key]
                enum += 1
        if enum > 2:
            node, elem, eroi = splitelem(node, elem, eroi, e2r, i)
            i += 3  # replace the old element with four new elements
        i += 1
    return node, elem, eroi


def demo_immc_vessel(cfg=None):
    if cfg is None:
        vs = loadmat(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../mmclab/example/vessel.mat"))
        edge = vs["edge"].astype(int)  # vessel edges
        fbox = vs["fbox"].astype(int)  # bounding box faces
        nbox = vs["nbox"].astype(float)  # bounding box nodes
        node = vs["node"].astype(float)  # nodes on the vessel network
        noder = np.ravel(vs["noder"]).astype(float)  # nodal radius of each node

        offset = nbox.shape[0]
        edge = edge + offset  # change the node index after inserting bounding box
        noder = np.concatenate([np.zeros(offset), noder])  # add offset into noder
        node = np.vstack([nbox, node])  # add nodes of bounding box into node
        fedge = np.column_stack([edge[:, 0], edge])  # convert edge to face for mesh generation
        face = [list(f) for f in np.vstack([fbox, fedge])]  # add faces of bounding box into face

        vnode, velem = i2m.s2m(node, face, 1, 100, "tetgen1.5", None, None, "-YY")[:2]
        velem = velem[:, :4].astype(int)
        e2r = edge2radius(edge, noder)
        vnode, velem, eroi = edgeroi(vnode, velem, e2r)

        res = 0.001  # resolution: 0.001mm
        cfg = {}
        cfg["nphoton"] = 1e6
        cfg["srctype"] = "pencil"
        cfg["srcpos"] = list(np.array([150.5, 150.5, 605]) * res)
        cfg["srcdir"] = [0, 0, -1]
        cfg["elemprop"] = np.ones(velem.shape[0])
        cfg["prop"] = [[0, 0, 1, 1], [0.0458, 35.6541, 0.9000, 1.3700], [23.0543, 9.3985, 0.9000, 1.3700]]  # dermis, blood
        cfg["tstart"] = 0
        cfg["tend"] = 5e-9
        cfg["tstep"] = 5e-9
        cfg["debuglevel"] = "TP"
        cfg["method"] = "grid"
        cfg["steps"] = [0.001, 0.001, 0.001]
        cfg["gpuid"] = -1
        cfg["isreflect"] = 1
        cfg["outputtype"] = "energy"
        cfg["elem"] = velem
        cfg["node"] = vnode * res
        cfg["edgeroi"] = eroi * res
        cfg["noderoi"] = np.zeros(vnode.shape[0])
        cfg["noderoi"][: len(noder)] = noder * res

    tic = time.time()
    flux = pmmc.mmclab(cfg)
    print("elapsed %.2f s" % (time.time() - tic))
    return cfg, flux


if __name__ == "__main__":
    cfg, flux = demo_immc_vessel()
    plt.imshow(np.log10(np.rot90(np.squeeze(flux["flux"][:, 149, :]))))
    plt.colorbar()
    plt.title("log10 energy deposit at y=150")
    plt.show()
