#!/usr/bin/env python3
"""
mkmeshcase.py - convert an MCX box-domain trajectory test case to an MMC input

Usage:
    mkmeshcase.py MCXCASE.json MMCCASE.json

MCXCASE uses a Grid shape covering Domain.Dim, optionally split along z by
ZLayers. The same box is tessellated into tetrahedra (6 Kuhn tetrahedra per
box of a rectilinear grid, conforming across boxes) and written inline as
Mesh.MeshNode/MeshElem, so that MMC traces photons in exactly the same domain
and media as MCX.

The box is also split along x and y at unequal fractions (XYSPLIT): with a
single box, the symmetric pencil beam at the domain center would travel inside
the internal x=y face shared by the Kuhn tetrahedra, a degenerate ray that
the tracer handles by nudging the photon (fixphoton), differently per backend.
"""

import json
import sys
from itertools import permutations


XYSPLIT = (0.37, 0.61)  # x and y split fractions, unequal to avoid diagonal faces through the center


def kuhn_tets(x0, x1, y0, y1, z0, z1):
    """6 tetrahedra filling [x0,x1]x[y0,y1]x[z0,z1], all sharing the main diagonal"""
    nodes, tets = [], []
    corner = {}

    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                corner[(i, j, k)] = len(nodes)
                nodes.append([x0 + i * (x1 - x0), y0 + j * (y1 - y0), z0 + k * (z1 - z0)])

    for perm in permutations(range(3)):
        v = [0, 0, 0]
        tet = [corner[tuple(v)]]

        for ax in perm:
            v[ax] = 1
            tet.append(corner[tuple(v)])

        tets.append(tet)

    return nodes, tets


def find_elem(nodes, elems, pos, direction):
    """1-based index of a tetrahedron enclosing pos nudged along direction"""
    import numpy as np

    p = np.asarray(pos[:3], dtype=float) + 1e-4 * np.asarray(direction[:3], dtype=float)

    for i, e in enumerate(elems):
        v = np.asarray([nodes[k - 1] for k in e[:4]], dtype=float)
        bary = np.linalg.solve(np.vstack([v.T, np.ones(4)]), np.append(p, 1.0))

        if bary.min() >= -1e-9:
            return i + 1

    raise ValueError("source is outside of the mesh")


def main():
    case = json.load(open(sys.argv[1]))
    dim = case["Domain"]["Dim"]
    tag = 1
    layers = [[0, dim[2], 1]]

    for shape in case.get("Shapes", []):
        if "Grid" in shape:
            tag = shape["Grid"].get("Tag", 1)
            layers = [[0, dim[2], tag]]
        elif "ZLayers" in shape:
            for zl in shape["ZLayers"]:  # [first, last, tag] in 1-based voxel index
                z0, z1, t = zl[0] - 1, zl[1], zl[2]
                # split the existing layer that contains [z0, z1]
                newlayers = []

                for a, b, lt in layers:
                    if b <= z0 or a >= z1:
                        newlayers.append([a, b, lt])
                        continue

                    if a < z0:
                        newlayers.append([a, z0, lt])

                    newlayers.append([max(a, z0), min(b, z1), t])

                    if b > z1:
                        newlayers.append([z1, b, lt])

                layers = newlayers

    nodemap, nodes, elems = {}, [], []

    xs = [0, dim[0] * XYSPLIT[0], dim[0]]
    ys = [0, dim[1] * XYSPLIT[1], dim[1]]
    boxes = [(xs[i], xs[i + 1], ys[j], ys[j + 1], z0, z1, t)
             for z0, z1, t in layers for i in range(2) for j in range(2)]

    for x0, x1, y0, y1, z0, z1, t in boxes:
        n, tets = kuhn_tets(x0, x1, y0, y1, z0, z1)
        local = []

        for p in n:
            key = tuple(p)

            if key not in nodemap:
                nodemap[key] = len(nodes)
                nodes.append(list(p))

            local.append(nodemap[key])

        for tet in tets:
            elems.append([local[i] + 1 for i in tet] + [t])

    optode = case["Optode"]
    e0 = find_elem(nodes, elems, optode["Source"]["Pos"], optode["Source"]["Dir"])

    out = {
        "Session": dict(case["Session"]),
        "Forward": case["Forward"],
        "Domain": {"Media": case["Domain"]["Media"], "LengthUnit": case["Domain"].get("LengthUnit", 1)},
        "Mesh": {"MeshID": case["Session"]["ID"], "InitElem": e0, "MeshNode": nodes, "MeshElem": elems},
        "Optode": optode,
    }
    json.dump(out, open(sys.argv[2], "w"), indent=1)


if __name__ == "__main__":
    main()
