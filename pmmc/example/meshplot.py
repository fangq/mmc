"""
Plotting helper shared by the pmmc examples: show a nodal quantity on a plane cut through a
tetrahedral mesh, the counterpart of plotmesh(..., 'y=30') and qmeshcut in the mmclab demos.
"""

import numpy as np
import matplotlib.pyplot as plt
import iso2mesh as i2m

_AXES = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}


def plotcut(node, elem, value, plane, ax=None, title=None, levels=40, **kwargs):
    """
    plotcut(node, elem, value, plane, ax=None, title=None)

    node, elem: the mesh (elem is 1-based, only its first 4 columns are used)
    value: a nodal value, length of node (or longer, the nodes of the mesh come first)
    plane: a plane such as 'y=30'
    """
    node = np.asarray(node)[:, :3]
    value = np.ravel(np.asarray(value))[: node.shape[0]]
    cutpos, cutval = i2m.qmeshcut(np.asarray(elem)[:, :4], node, value, plane)[:2]
    cutpos, cutval = np.asarray(cutpos), np.ravel(cutval)
    good = np.isfinite(cutval)
    i, j = _AXES[plane.strip()[0]]
    ax = ax or plt.gca()
    tc = ax.tricontourf(cutpos[good, i], cutpos[good, j], cutval[good], levels, **kwargs)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title, fontsize=9)
    plt.colorbar(tc, ax=ax)
    return tc
