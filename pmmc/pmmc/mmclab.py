# Copyright (c) 2026 Qianqian Fang <q.fang at neu.edu>. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""pmmc.mmclab() and the mesh helpers it uses, ported from mmclab/mmclab.m and matlab/mmc*.m

mmclab() prepares the mesh and the source as mmclab.m does (element reorientation, face
neighbors, element volumes, the initial element, and the source domain of wide-field sources),
runs pmmc.run(), and unpacks the outputs the same way. The mesh functions need iso2mesh
(pip install iso2mesh).
"""

import copy
import numpy as np

__all__ = ["mmclab", "mmcsrcdomain", "mmcaddsrc", "mmcadddet", "mmcraytrace"]

WIDEFIELD_SRC = ("planar", "pattern", "fourier", "fourierx", "fourierx2d")


def _i2m():
    try:
        import iso2mesh
    except ImportError:  # pragma: no cover
        raise ImportError('pmmc.mmclab needs iso2mesh to prepare the mesh, please run "pip install iso2mesh"')
    return iso2mesh


def _first(x):
    """the 1st output of an iso2mesh function, which may return a tuple"""
    return x[0] if isinstance(x, tuple) else x


def _run(cfg):
    from _pmmc import run

    return run(cfg)


def mmcsrcdomain(cfg, meshbbx, expansion=1.1, rotate=0.0, keepshape=0, circlediv=100):
    """
    srcnode, srcface = mmcsrcdomain(cfg, meshbbx, expansion=1.1, rotate=0.0, keepshape=0)

    Define a source domain (a triangle or a polygon on the source plane) for wide-field sources
    outside of the mesh, ported from mmcsrcdomain.m

    input:
        cfg: a dict with at least srctype, srcpos and srcdir, and srcparam1/srcparam2 as needed
        meshbbx: the bounding box of the mesh, [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        expansion: the expansion ratio of the source domain, 1.1 by default
        rotate: the rotation angle (in radian) of the domain around srcdir
        keepshape: 1 to keep the source shape (planar/disk) instead of a covering triangle
        circlediv: the number of divisions of a disk source domain when keepshape is 1

    output:
        srcnode: the nodes of the source domain, Nx3
        srcface: the triangles of the source domain, 1-based
    """
    i2m = _i2m()
    for key in ("srcpos", "srcdir", "srctype"):
        if key not in cfg:
            raise ValueError("you must at least define cfg['srcpos'], cfg['srcdir'], cfg['srctype']")

    srctype = cfg["srctype"]
    srcpos = np.ravel(np.asarray(cfg["srcpos"], dtype=float))[:3]
    srcdir = np.ravel(np.asarray(cfg["srcdir"], dtype=float)).copy()
    meshbbx = np.asarray(meshbbx, dtype=float)

    if srctype == "isotropic":
        srcdir[:3] = np.random.rand(3)
    srcdir[:3] = srcdir[:3] / np.linalg.norm(srcdir[:3])

    domainradius = np.min(meshbbx[1] - meshbbx[0]) * 0.25
    domaindiagnoal = np.linalg.norm(meshbbx[1] - meshbbx[0])
    vperp = np.cross(srcdir[:3], np.random.rand(3))
    param1 = np.ravel(np.asarray(cfg.get("srcparam1", [0, 0, 0, 0]), dtype=float))
    param2 = np.ravel(np.asarray(cfg.get("srcparam2", [0, 0, 0, 0]), dtype=float))

    if srctype in ("pencil", "cone", "zgaussian", "isotropic", "arcsine"):
        srcnode = i2m.orthdisk(srcpos, srcpos + srcdir[:3], domainradius, 3, vperp, rotate)
        srcface = np.array([[1, 2, 3]])
    elif srctype == "gaussian":
        srcnode = i2m.orthdisk(srcpos, srcpos + srcdir[:3], domaindiagnoal * 0.5, 3, vperp, rotate)
        srcface = np.array([[1, 2, 3]])
    elif srctype in ("line", "slit"):
        v0 = srcpos + param1[:3] * 0.5
        srcnode = i2m.orthdisk(v0, v0 + srcdir[:3], np.linalg.norm(param1[:3]) * expansion, 3, param1[:3])
        srcface = np.array([[1, 2, 3]])
    elif srctype in WIDEFIELD_SRC:
        v0, v1, v2 = srcpos, param1[:3], param2[:3]
        if srctype in ("fourierx", "fourierx2d"):
            v2 = np.cross(srcdir[:3], v1) * (param1[3] / np.linalg.norm(v1))
        voff1 = (v1 + v2) * 0.5
        voff2 = (v1 - v2) * 0.5
        if keepshape:
            srcnode = np.array([v0, v0 + v1, v0 + v1 + v2, v0 + v2]) + np.array([-voff1, voff2, voff1, -voff2]) * (expansion - 1.0)
            srcface = np.array([[1, 2, 3], [3, 4, 1]])
        else:
            radius = max(np.linalg.norm(voff1), np.linalg.norm(voff2)) * 2.0 * expansion
            srcnode = i2m.orthdisk(srcpos + voff1, srcpos + voff1 + srcdir[:3], radius, 3, v1 + v2, rotate)
            srcface = np.array([[1, 2, 3]])
    elif srctype == "disk":
        if keepshape:
            from scipy.spatial import Delaunay

            srcnode = i2m.orthdisk(srcpos, srcpos + srcdir[:3], param1[0] * expansion, circlediv)
            srcface = Delaunay(srcnode[:, :2]).simplices + 1
        else:
            srcnode = i2m.orthdisk(srcpos, srcpos + srcdir[:3], param1[0] * 2.0 * expansion, 3, vperp, rotate)
            srcface = np.array([[1, 2, 3]])
    else:
        raise ValueError("source type not supported: " + str(srctype))

    return np.asarray(srcnode, dtype=float), srcface


def mmcaddsrc(node, elem, src, **opt):
    """
    newnode, newelem = mmcaddsrc(node, elem, src, **opt)

    Add a wide-field source domain to the mesh, ported from mmcaddsrc.m. The space between the
    mesh and the source domain is tessellated and labeled 0, and the elements attached to the
    source domain are labeled -1 (extcorelabel)

    input:
        node, elem: the mesh, elem is 1-based, with an optional 5th column of labels
        src: the nodes of the source domain (Nx3), or a source dict passed to mmcsrcdomain
        opt: extcmdopt (tetgen options, '-Y' by default), extlabel, extcorelabel, and the
             mmcsrcdomain options (expansion, rotate, keepshape) when src is a dict

    output:
        newnode, newelem: the merged mesh, newelem has 5 columns
    """
    i2m = _i2m()
    node = np.asarray(node, dtype=float)
    elem = np.asarray(elem)
    opt = {k.lower(): v for k, v in opt.items()}
    opt.setdefault("extcmdopt", "-Y")
    if isinstance(src, dict):
        domopt = {k: opt[k] for k in ("expansion", "rotate", "keepshape", "circlediv") if k in opt}
        src = mmcsrcdomain(src, [node[:, :3].min(axis=0), node[:, :3].max(axis=0)], **domopt)[0]
    opt["newnode"] = np.asarray(src, dtype=float)
    if elem.shape[1] == 4:
        elem = np.hstack([elem, np.ones((elem.shape[0], 1), dtype=elem.dtype)])
    # iso2mesh-python reads newnode/extcmdopt from keyword arguments and extlabel/extcorelabel
    # from the option dict
    newnode, newelem = i2m.meshrefine(node, elem, opt, **opt)[:2]
    return newnode, newelem


def mmcadddet(node, elem, det, **opt):
    """
    newnode, newelem = mmcadddet(node, elem, det, **opt)

    Add a wide-field detector domain to the mesh, labeled -2, ported from mmcadddet.m
    """
    opt = dict(opt, extcorelabel=-2, keepshape=1, expansion=1.0)
    return mmcaddsrc(node, elem, det, **opt)


def mmcraytrace(node, elem, p0, v0, e0="-"):
    """
    p, e0 = mmcraytrace(node, elem, p0, v0, e0)

    Find the entry point of a ray on the mesh surface and the element enclosing it, ported from
    mmcraytrace.m

    input:
        node, elem: the mesh (tetrahedra or a triangular surface), 1-based
        p0, v0: the origin and the direction of the ray
        e0: '>' searches along v0, '<' against v0, '-' or nan in both directions

    output:
        p: the intersection closest to p0 in the searched direction
        e0: the (1-based) element or surface triangle containing p
    """
    i2m = _i2m()
    node = np.asarray(node, dtype=float)
    elem = np.asarray(elem)
    p0 = np.ravel(np.asarray(p0, dtype=float))[:3]
    v0 = np.ravel(np.asarray(v0, dtype=float))[:3]
    face = elem[:, :3] if elem.shape[1] == 3 else _first(i2m.volface(elem[:, :4]))
    t, u, v, idx = i2m.raytrace(p0, v0, node[:, :3], face)
    if len(idx) == 0:
        raise ValueError("ray does not intersect with the mesh")
    t = np.asarray(t)[np.asarray(idx) - 1]  # raytrace returns 1-based triangle indices
    if isinstance(e0, str) and e0 == ">":
        idx1 = np.where(t >= 0)[0]
    elif isinstance(e0, str) and e0 == "<":
        idx1 = np.where(t <= 0)[0]
    elif (isinstance(e0, str) and e0 == "-") or (not isinstance(e0, str) and np.isnan(e0)):
        idx1 = np.arange(len(t))
    else:
        raise ValueError("ray direction specifier is not recognized")
    if len(idx1) == 0:
        raise ValueError("no intersection is found along the ray direction")
    loc = np.argmin(np.abs(t[idx1]))
    faceidx = idx[idx1[loc]]
    p = p0 + t[idx1[loc]] * v0
    if elem.shape[1] == 3:
        return p, int(faceidx)
    # the element that has this triangle as one of its faces
    felem = np.sort(face[faceidx - 1, :3])
    tet = elem[:, :4]
    allface = np.vstack([tet[:, [0, 1, 2]], tet[:, [1, 0, 3]], tet[:, [0, 2, 3]], tet[:, [1, 3, 2]]])
    match = np.where(np.all(np.sort(allface, axis=1) == felem, axis=1))[0]
    loc = (match[0] % tet.shape[0]) + 1
    return p, int(loc)


def _reorient(cfg, i2m):
    """reorient the elements to positive volume, permuting the edge/face ROIs of flipped elements"""
    elem, _, idx = i2m.meshreorient(cfg["node"][:, :3], np.array(cfg["elem"][:, :4]))
    cfg["elem"] = elem
    if "edgeroi" in cfg:
        cfg["edgeroi"] = np.array(cfg["edgeroi"])
        cfg["edgeroi"][idx, :] = cfg["edgeroi"][np.ix_(idx, [0, 2, 1, 4, 3, 5])]
    if "faceroi" in cfg:
        cfg["faceroi"] = np.array(cfg["faceroi"])
        cfg["faceroi"][idx, :] = cfg["faceroi"][np.ix_(idx, [1, 0, 2, 3])]
    cfg["isreoriented"] = 1


def _tsearch(node, elem, pos):
    """the 1-based enclosing element of a point, nan if outside, as tsearchn in MATLAB/Octave"""
    idx = _first(_i2m().tsearchn(node[:, :3], elem[:, :4], np.atleast_2d(pos)[:, :3]))
    return float(np.ravel(idx)[0])


def _prepare(cfg):
    """the mesh and source preparation of mmclab.m for one simulation"""
    i2m = _i2m()
    cfg = copy.copy(cfg)
    if "node" not in cfg or "elem" not in cfg:
        raise ValueError("cfg['node'] or cfg['elem'] is missing")
    for key in ("srcpattern", "srcpos", "detpos", "prop", "workload", "srcdir"):
        if key in cfg:
            cfg[key] = np.array(cfg[key], dtype=np.float64)
    cfg["node"] = np.asarray(cfg["node"], dtype=np.float64)
    cfg["elem"] = np.asarray(cfg["elem"])
    if cfg["elem"].shape[1] > 4:
        cfg["elemprop"] = cfg["elem"][:, 4]
    if not cfg.get("isreoriented"):
        _reorient(cfg, i2m)
    if cfg.get("facenb") is None or np.size(cfg["facenb"]) == 0:
        cfg["facenb"] = _first(i2m.faceneighbors(cfg["elem"][:, :4]))
    if cfg.get("evol") is None or np.size(cfg["evol"]) == 0:
        cfg["evol"] = i2m.elemvolume(cfg["node"][:, :3], cfg["elem"][:, :4])
    isoptix = cfg.get("compute") == "optix"
    if not isoptix and np.any(np.asarray(cfg["evol"]) == 0):
        raise ValueError(
            "input mesh can not contain degenerated elements (%s), please double check your input mesh; "
            "if you use a widefield source, please rerun mmcsrcdomain and set expansion to a larger value"
            % (np.where(np.asarray(cfg["evol"]) == 0)[0] + 1)
        )
    if "srcpos" not in cfg:
        raise ValueError("cfg['srcpos'] field is missing")
    if "srcdir" not in cfg:
        raise ValueError("cfg['srcdir'] field is missing")
    srcpos = np.atleast_2d(cfg["srcpos"])
    srcdir = np.atleast_2d(cfg["srcdir"])
    if cfg.get("e0") is None:
        # only the first row of a multi-source srcpos is looked up, mmc finds the others
        cfg["e0"] = _tsearch(cfg["node"], cfg["elem"], srcpos[0, :3])
    e0 = cfg["e0"]
    if (not isinstance(e0, str) and np.isnan(e0) and cfg.get("srctype") == "pencil") or isinstance(e0, str):
        print("searching initial element ...")
        srcpos = np.array(srcpos)
        srcpos[0, :3], cfg["e0"] = mmcraytrace(cfg["node"], cfg["elem"], srcpos[0, :3], srcdir[0, :3], e0)
        cfg["srcpos"] = srcpos if np.ndim(cfg["srcpos"]) > 1 else srcpos[0]
    # a 2-D srcpattern is one pattern; a 3-D one is a [srcnum, Nx, Ny] stack (the reshape of a 2-D
    # pattern in mmclab.m to Nx x Ny x 1 is a no-op in MATLAB and must not add an axis here)

    # wide-field sources pick a launch position on the source plane for every photon, so the
    # launch element is searched in the source domain (elements labeled -1) added to the mesh;
    # this is keyed off the source type, as a srcpos on the mesh surface still has a valid e0
    # (see mmclab.m)
    iswidesrc = cfg.get("srctype") in WIDEFIELD_SRC and not isoptix
    if (not isinstance(cfg["e0"], str) and np.isnan(cfg["e0"])) or iswidesrc:
        if "srcparam1" not in cfg or "srcparam2" not in cfg:
            raise ValueError("for wide-field sources, you must provide srcparam1 and srcparam2")
        cfg.setdefault("srctype", "pencil")
        srcdef = {k: cfg[k] for k in ("srctype", "srcpos", "srcdir", "srcparam1", "srcparam2")}
        node = cfg["node"]
        sdom = mmcsrcdomain(srcdef, [node[:, :3].min(axis=0), node[:, :3].max(axis=0)])[0]
        rounded = {tuple(r) for r in np.round(node[:, :3] * 1e10) * 1e-10}
        isinside = [tuple(r) in rounded for r in np.round(sdom * 1e10) * 1e-10]
        if not any(isinside) and not isoptix:
            elem = cfg["elem"]
            if elem.shape[1] == 4:
                elem = np.hstack([elem, np.reshape(cfg["elemprop"], (-1, 1)).astype(elem.dtype)])
            cfg["node"], elem = mmcaddsrc(cfg["node"], elem, sdom)
            cfg["elemprop"] = elem[:, 4]
            cfg["elem"] = elem
            _reorient(cfg, i2m)
            cfg["facenb"] = _first(i2m.faceneighbors(cfg["elem"][:, :4]))
            cfg["evol"] = i2m.elemvolume(cfg["node"][:, :3], cfg["elem"][:, :4])
        elemprop = np.ravel(cfg["elemprop"])
        if cfg["srctype"] in ("pencil", "isotropic", "cone", "zgaussian"):
            cfg["e0"] = _tsearch(cfg["node"], cfg["elem"], cfg["srcpos"])
            if np.isnan(cfg["e0"]):
                cfg["e0"] = -1
        elif np.count_nonzero(elemprop < 0) == 1:
            cfg["e0"] = int(np.where(elemprop < 0)[0][0]) + 1
        else:
            cfg["e0"] = -1
        if iswidesrc and not np.any(elemprop < 0):
            raise ValueError(
                "no source-domain element was tessellated into the mesh for the '%s' wide-field source, "
                "so the per-photon launch-element search has nothing to search and the simulation would "
                "deposit the entire beam into element cfg['e0']. This normally means cfg['srcpos'] sits on "
                "or inside the mesh surface: move the source outside the mesh, against cfg['srcdir'], so "
                "that mmcaddsrc can attach the source domain" % cfg["srctype"]
            )
    if "elemprop" not in cfg:
        raise ValueError("cfg['elemprop'] field is missing")
    if "nphoton" not in cfg:
        raise ValueError("cfg['nphoton'] field is missing")
    elemprop = np.ravel(cfg["elemprop"])
    if "prop" not in cfg or np.shape(cfg["prop"])[0] < np.max(elemprop) + 1 or np.all(elemprop <= 0):
        raise ValueError("cfg['prop'] field is missing or insufficient")
    if not isinstance(cfg["e0"], str):
        cfg["e0"] = int(cfg["e0"])
    return cfg


def _unpackdetp(detp, cfg):
    """split the detected photon records into named fields, as mmclab.m"""
    medianum = np.shape(cfg["prop"])[0] - 1
    newdetp = {}
    newdetp["detid"] = detp[0, :].astype(np.int32)
    newdetp["nscat"] = detp[1 : medianum + 1, :].astype(np.int32).T  # number of scattering events per medium
    newdetp["ppath"] = detp[medianum + 1 : 2 * medianum + 1, :].T  # partial path lengths per medium
    if cfg.get("ismomentum"):
        newdetp["mom"] = detp[2 * medianum + 1 : 3 * medianum + 1, :].T  # momentum transfer per medium
    if cfg.get("issaveexit"):
        newdetp["p"] = detp[-7:-4, :].T  # exit positions
        newdetp["v"] = detp[-4:-1, :].T  # exit directions
    newdetp["w0"] = detp[-1, :].copy()  # the initial packet weight
    if np.ndim(cfg.get("srcpattern")) == 3:
        newdetp["w0"] = newdetp["w0"].astype(np.float32).view(np.uint32)
    newdetp["prop"] = cfg["prop"]
    if "unitinmm" in cfg:
        newdetp["unitinmm"] = cfg["unitinmm"]
    newdetp["data"] = detp
    return newdetp


def _unpacktraj(data):
    """sort the trajectory records by photon id and split them into named fields, as mmclab.m"""
    data = np.asarray(data, dtype=np.float32)
    pid = data[0, :].view(np.uint32)
    idx = np.argsort(pid, kind="stable")
    traj = {}
    traj["id"] = pid[idx]
    traj["pos"] = data[1:4, idx].T
    traj["eid"] = data[5, idx].view(np.uint32)
    traj["data"] = np.vstack([traj["id"].astype(np.float32), data[1:, idx]])
    return traj


def mmclab(cfg, option=None):
    """
    res = pmmc.mmclab(cfg)
    res = pmmc.mmclab(cfg, 'prep')
    gpuinfo = pmmc.mmclab('gpuinfo')

    Run a mesh-based Monte Carlo (MMC) simulation with the same pre- and post-processing as
    mmclab.m; pmmc.run() only runs the simulation

    input:
        cfg: a dict of simulation settings, or a list of such dicts, see the mmclab documentation;
             node/elem/elemprop (or elem with a 5th label column), prop, srcpos, srcdir, nphoton,
             tstart/tend/tstep are needed. Before the simulation, mmclab
               - reorients elements to a positive volume (unless cfg['isreoriented'] is 1),
                 and computes cfg['facenb'] and cfg['evol'] if missing
               - finds the initial element cfg['e0'] (raytracing to the mesh surface for pencil
                 beams outside the mesh, or when e0 is '>', '<' or '-')
               - adds a source domain to the mesh for wide-field sources (planar, pattern,
                 fourier, fourierx, fourierx2d) outside the mesh, see mmcsrcdomain/mmcaddsrc
             mmclab.m sets the saved outputs from its number of outputs; here they are set by
             cfg['issavedet'] (detected photons), cfg['issaveseed'] (their seeds) and
             cfg['debuglevel'] containing 'M' (trajectories)
        option: 'prep' returns the prepared cfg without simulating

    output:
        res: a dict (or a list of dicts for a list of cfg) with the outputs of pmmc.run, where
             res['flux']: the fluence rate (or the output selected by cfg['outputtype']) at the
                    nodes or elements, a complex array for RF (omega>0) simulations
             res['dref']: diffuse reflectance on the surface faces, if cfg['issaveref'] is 1
             res['jmua'], res['jd'], ...: adjoint Jacobians for the adjoint output types
             res['detp']: detected photons, a dict with detid, nscat, ppath, mom, p, v, w0,
                    prop, unitinmm and data (the raw records), as mmclab.m returns
             res['seeds']: the RNG seeds of the detected photons, one column per photon
             res['traj']: photon trajectories sorted by photon id, a dict with id, pos, eid, data
             res['cfg']: the prepared cfg, whose node/elem include any added source domain;
                    nodal outputs follow this mesh (the 3rd output of mmclab.m)
    """
    if isinstance(cfg, str):
        from _pmmc import gpuinfo, version

        if cfg == "gpuinfo":
            return gpuinfo()
        if cfg == "version":
            return version()
        raise ValueError("unsupported command: " + cfg)
    if isinstance(cfg, (list, tuple)):
        return [mmclab(c, option) for c in cfg]
    if not isinstance(cfg, dict):
        raise ValueError("cfg must be a dict or a list of dicts")

    cfg = _prepare(cfg)
    if option == "prep":
        return cfg
    if option not in (None, "omp", "sse", "opencl", "cuda"):
        raise ValueError("option is not recognized: " + str(option))

    if isinstance(cfg.get("gpuid"), (int, np.integer)) and cfg["gpuid"] < -1:
        cfg["gpuid"] = 1
    res = _run(cfg)

    if "fluximag" in res:  # RF output: return the fluence as a complex array, as mmclab.m does
        res["flux"] = res["flux"] + 1j * np.asarray(res.pop("fluximag"), dtype=np.float64)
    for key in [k for k in res if k.endswith("_im") and k[:-3] in res]:  # RF adjoint Jacobians (jmua_im, jd_im, ...)
        res[key[:-3]] = res[key[:-3]] + 1j * np.asarray(res.pop(key))
    if "detp" in res and cfg.get("issaveexit", 0) != 2 and np.size(res["detp"]) > 0:
        res["detp"] = _unpackdetp(np.asarray(res["detp"]), cfg)
    if "traj" in res and np.size(res["traj"]) > 0:
        res["traj"] = _unpacktraj(res["traj"])
    res["cfg"] = cfg
    return res
