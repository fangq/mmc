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

"""The mesh front end of pmmc: :func:`mmclab` and its helpers, ported from ``mmclab.m``
and ``mmcsrcdomain.m``, ``mmcaddsrc.m``, ``mmcadddet.m`` and ``mmcraytrace.m``.

:func:`mmclab` prepares the mesh and the source as ``mmclab.m`` does (element
reorientation, face neighbors, element volumes, the initial element, and the source domain
of wide-field sources), runs ``pmmc.run``, and unpacks the outputs the same way. The mesh
functions need ``iso2mesh`` (``pip install iso2mesh``), which is imported only when they
are called.

Element and node indices in the inputs and outputs are 1-based, as in MATLAB.
"""

import copy
import numpy as np

__all__ = ["mmclab", "mmcsrcdomain", "mmcaddsrc", "mmcadddet", "mmcraytrace"]

WIDEFIELD_SRC = ("planar", "pattern", "fourier", "fourierx", "fourierx2d")


def _i2m():
    """Import and return the ``iso2mesh`` module, with an install hint if it is missing."""
    try:
        import iso2mesh
    except ImportError:  # pragma: no cover
        raise ImportError(
            'pmmc.mmclab needs iso2mesh to prepare the mesh, please run "pip install iso2mesh"'
        )
    return iso2mesh


def _first(x):
    """the 1st output of an iso2mesh function, which may return a tuple"""
    return x[0] if isinstance(x, tuple) else x


def _run(cfg):
    """Run ``_pmmc.run(cfg)``; imported here so that the mesh helpers work without ``_pmmc``."""
    from _pmmc import run

    return run(cfg)


def mmcsrcdomain(cfg, meshbbx, expansion=1.1, rotate=0.0, keepshape=0, circlediv=100):
    """Define the source domain of a source outside the mesh, ported from ``mmcsrcdomain.m``.

    The domain is a triangle (or, with ``keepshape``, the source's own polygon) on the source
    plane, perpendicular to ``srcdir``, that covers the source aperture. :func:`mmcaddsrc`
    tessellates it into the mesh, so that photons launched anywhere on the aperture start in
    a known element.

    Parameters
    ----------
    cfg : dict
        The source: ``srctype``, ``srcpos`` and ``srcdir``, and ``srcparam1``/``srcparam2``
        as the source type needs. Supported types: pencil, cone, zgaussian, isotropic,
        arcsine, gaussian, line, slit, planar, pattern, fourier, fourierx, fourierx2d and
        disk.
    meshbbx : array_like, shape (2, 3)
        The bounding box of the mesh, ``[[xmin, ymin, zmin], [xmax, ymax, zmax]]``.
    expansion : float, optional
        The ratio by which the domain is enlarged beyond the aperture (default 1.1).
    rotate : float, optional
        The rotation of the triangle around ``srcdir``, in radians (default 0).
    keepshape : int, optional
        1 to use the shape of the aperture (a parallelogram for planar/pattern/fourier
        sources, a polygon for a disk) instead of a covering triangle (default 0).
    circlediv : int, optional
        The number of segments of a disk domain with ``keepshape`` (default 100).

    Returns
    -------
    srcnode : ndarray, shape (n, 3)
        The nodes of the source domain.
    srcface : ndarray, shape (m, 3)
        Its triangles, 1-based.

    Raises
    ------
    ValueError
        If ``srcpos``, ``srcdir`` or ``srctype`` is missing, or the type is not supported.
    """
    i2m = _i2m()
    for key in ("srcpos", "srcdir", "srctype"):
        if key not in cfg:
            raise ValueError(
                "you must at least define cfg['srcpos'], cfg['srcdir'], cfg['srctype']"
            )

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
        srcnode = i2m.orthdisk(
            srcpos, srcpos + srcdir[:3], domainradius, 3, vperp, rotate
        )
        srcface = np.array([[1, 2, 3]])
    elif srctype in ("gaussian", "hyperboloid"):
        # the beam radius at the source plane: the waist param1[0], widened for a hyperboloid
        # param1[1] from its focus (Rayleigh range param1[2]); the triangle, of circumradius at
        # least 8 beam radii, holds 4 of them, so the Gaussian tails, which launched photons
        # outside of the domain (a fatal error), have a chance of about 1e-14
        beamradius = param1[0]
        if srctype == "hyperboloid" and param1[2] > 0:
            beamradius *= np.sqrt(1 + (param1[1] / param1[2]) ** 2)
        srcnode = i2m.orthdisk(
            srcpos,
            srcpos + srcdir[:3],
            max(domaindiagnoal, 8 * beamradius),
            3,
            vperp,
            rotate,
        )
        srcface = np.array([[1, 2, 3]])
    elif srctype in ("line", "slit"):
        v0 = srcpos + param1[:3] * 0.5
        srcnode = i2m.orthdisk(
            v0, v0 + srcdir[:3], np.linalg.norm(param1[:3]) * expansion, 3, param1[:3]
        )
        srcface = np.array([[1, 2, 3]])
    elif srctype in WIDEFIELD_SRC or srctype == "pencilarray":
        v0, v1, v2 = srcpos, param1[:3], param2[:3]
        if srctype in ("fourierx", "fourierx2d"):
            v2 = np.cross(srcdir[:3], v1) * (param1[3] / np.linalg.norm(v1))
        voff1 = (v1 + v2) * 0.5
        voff2 = (v1 - v2) * 0.5
        if keepshape:
            srcnode = np.array([v0, v0 + v1, v0 + v1 + v2, v0 + v2]) + np.array(
                [-voff1, voff2, voff1, -voff2]
            ) * (expansion - 1.0)
            srcface = np.array([[1, 2, 3], [3, 4, 1]])
        else:
            radius = max(np.linalg.norm(voff1), np.linalg.norm(voff2)) * 2.0 * expansion
            srcnode = i2m.orthdisk(
                srcpos + voff1, srcpos + voff1 + srcdir[:3], radius, 3, v1 + v2, rotate
            )
            srcface = np.array([[1, 2, 3]])
    elif srctype in ("disk", "ring"):
        radius = param1[0]  # the outer radius
        if keepshape:
            from scipy.spatial import Delaunay

            srcnode = i2m.orthdisk(
                srcpos, srcpos + srcdir[:3], radius * expansion, circlediv
            )
            srcface = Delaunay(srcnode[:, :2]).simplices + 1
        else:
            srcnode = i2m.orthdisk(
                srcpos,
                srcpos + srcdir[:3],
                radius * 2.0 * expansion,
                3,
                vperp,
                rotate,
            )
            srcface = np.array([[1, 2, 3]])
    else:
        raise ValueError("source type not supported: " + str(srctype))

    return np.asarray(srcnode, dtype=float), srcface


def mmcaddsrc(node, elem, src, **opt):
    """Add a source domain to a mesh, ported from ``mmcaddsrc.m``.

    The space between the mesh and the source domain is tessellated with tetgen (through
    ``iso2mesh.meshrefine``) and labeled ``extlabel`` (0, background), and the elements
    touching the source domain are labeled ``extcorelabel`` (-1), where photons are launched.

    Parameters
    ----------
    node : array_like, shape (nn, 3)
        The mesh nodes.
    elem : array_like, shape (ne, 4) or (ne, 5)
        The tetrahedra, 1-based, with an optional 5th column of labels (1 if absent).
    src : array_like or dict
        The nodes of the source domain (``(n, 3)``), or a source dict, which is passed to
        :func:`mmcsrcdomain` with the bounding box of the mesh.
    **opt
        ``extcmdopt`` (tetgen options, default ``'-Y'``), ``extlabel``, ``extcorelabel``,
        and, when ``src`` is a dict, the ``expansion``, ``rotate``, ``keepshape`` and
        ``circlediv`` options of :func:`mmcsrcdomain`.

    Returns
    -------
    newnode : ndarray
        The nodes of the merged mesh.
    newelem : ndarray, shape (ne2, 5)
        Its tetrahedra and labels, 1-based.
    """
    i2m = _i2m()
    node = np.asarray(node, dtype=float)
    elem = np.asarray(elem)
    opt = {k.lower(): v for k, v in opt.items()}
    opt.setdefault("extcmdopt", "-Y")
    if isinstance(src, dict):
        domopt = {
            k: opt[k]
            for k in ("expansion", "rotate", "keepshape", "circlediv")
            if k in opt
        }
        src = mmcsrcdomain(
            src, [node[:, :3].min(axis=0), node[:, :3].max(axis=0)], **domopt
        )[0]
    opt["newnode"] = np.asarray(src, dtype=float)
    if elem.shape[1] == 4:
        elem = np.hstack([elem, np.ones((elem.shape[0], 1), dtype=elem.dtype)])
    # iso2mesh-python reads newnode/extcmdopt from keyword arguments and extlabel/extcorelabel
    # from the option dict
    newnode, newelem = i2m.meshrefine(node, elem, opt, **opt)[:2]
    return newnode, newelem


def mmcadddet(node, elem, det, **opt):
    """Add a wide-field detector domain to a mesh, ported from ``mmcadddet.m``.

    The same as :func:`mmcaddsrc` with ``extcorelabel=-2``, ``keepshape=1`` and
    ``expansion=1``: the elements touching the detector aperture are labeled -2.

    Parameters
    ----------
    node, elem : array_like
        The mesh, as in :func:`mmcaddsrc`.
    det : array_like or dict
        The nodes of the detector domain, or a dict with ``srctype``, ``srcpos``,
        ``srcdir``, ``srcparam1`` and ``srcparam2`` describing the detector aperture as a
        source.
    **opt
        Other options of :func:`mmcaddsrc`.

    Returns
    -------
    newnode, newelem : ndarray
        The merged mesh, as in :func:`mmcaddsrc`.
    """
    opt = dict(opt, extcorelabel=-2, keepshape=1, expansion=1.0)
    return mmcaddsrc(node, elem, det, **opt)


def mmcraytrace(node, elem, p0, v0, e0="-"):
    """Find where a ray meets the mesh surface, and the element there, ported from ``mmcraytrace.m``.

    Parameters
    ----------
    node : array_like, shape (nn, 3)
        The mesh nodes.
    elem : array_like
        The tetrahedra (``(ne, 4)`` or more columns), or a triangular surface (``(nf, 3)``),
        1-based.
    p0 : array_like
        The origin of the ray.
    v0 : array_like
        The direction of the ray.
    e0 : {'-', '>', '<'} or nan, optional
        Search along ``v0`` (``'>'``), against it (``'<'``), or both ways (``'-'`` or nan,
        the default).

    Returns
    -------
    p : ndarray, shape (3,)
        The intersection closest to ``p0`` in the searched direction.
    e0 : int
        The 1-based tetrahedron (or surface triangle) that contains ``p``.

    Raises
    ------
    ValueError
        If the ray misses the mesh, or ``e0`` is not recognized.
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
    elif (isinstance(e0, str) and e0 == "-") or (
        not isinstance(e0, str) and np.isnan(e0)
    ):
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
    allface = np.vstack(
        [tet[:, [0, 1, 2]], tet[:, [1, 0, 3]], tet[:, [0, 2, 3]], tet[:, [1, 3, 2]]]
    )
    match = np.where(np.all(np.sort(allface, axis=1) == felem, axis=1))[0]
    loc = (match[0] % tet.shape[0]) + 1
    return p, int(loc)


def _reorient(cfg, i2m):
    """Reorient the elements of ``cfg`` to positive volumes, permuting the ``edgeroi``/``faceroi``
    of the flipped elements, and set ``cfg['isreoriented']``.
    """
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
    """Return the 1-based element enclosing a point, or nan if it is outside, as ``tsearchn``."""
    idx = _first(_i2m().tsearchn(node[:, :3], elem[:, :4], np.atleast_2d(pos)[:, :3]))
    return float(np.ravel(idx)[0])


def _prepare(cfg):
    """Prepare the mesh and the source of one simulation as ``mmclab.m`` does; see :func:`mmclab`.
    Returns a shallow copy of ``cfg``.
    """
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
    if (
        not isinstance(e0, str) and np.isnan(e0) and cfg.get("srctype") == "pencil"
    ) or isinstance(e0, str):
        print("searching initial element ...")
        srcpos = np.array(srcpos)
        srcpos[0, :3], cfg["e0"] = mmcraytrace(
            cfg["node"], cfg["elem"], srcpos[0, :3], srcdir[0, :3], e0
        )
        cfg["srcpos"] = srcpos if np.ndim(cfg["srcpos"]) > 1 else srcpos[0]
    # a 2-D srcpattern is one pattern; a 3-D one is a [srcnum, Nx, Ny] stack (the reshape of a 2-D
    # pattern in mmclab.m to Nx x Ny x 1 is a no-op in MATLAB and must not add an axis here)

    # wide-field sources pick a launch position on the source plane for every photon, so the
    # launch element is searched in the source domain (elements labeled -1) added to the mesh;
    # this is keyed off the source type, as a srcpos on the mesh surface still has a valid e0
    # (see mmclab.m)
    iswidesrc = cfg.get("srctype") in WIDEFIELD_SRC and not isoptix
    if (not isinstance(cfg["e0"], str) and np.isnan(cfg["e0"])) or iswidesrc:
        # srcparam1 and srcparam2 span the area of the planar-type sources; the others (disk,
        # gaussian, line, ...) use only srcparam1, and those not using it default to zeros
        for key in ("srcparam1", "srcparam2"):
            if cfg.get(key) is None or np.size(cfg[key]) == 0:
                if iswidesrc:
                    raise ValueError(
                        "for wide-field sources, you must provide srcparam1 and srcparam2"
                    )
                cfg[key] = np.zeros(4)
        cfg.setdefault("srctype", "pencil")
        srcdef = {
            k: cfg[k] for k in ("srctype", "srcpos", "srcdir", "srcparam1", "srcparam2")
        }
        node = cfg["node"]
        sdom = mmcsrcdomain(srcdef, [node[:, :3].min(axis=0), node[:, :3].max(axis=0)])[
            0
        ]
        rounded = {tuple(r) for r in np.round(node[:, :3] * 1e10) * 1e-10}
        isinside = [tuple(r) in rounded for r in np.round(sdom * 1e10) * 1e-10]
        if not any(isinside) and not isoptix:
            elem = cfg["elem"]
            if elem.shape[1] == 4:
                elem = np.hstack(
                    [elem, np.reshape(cfg["elemprop"], (-1, 1)).astype(elem.dtype)]
                )
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
    if (
        "prop" not in cfg
        or np.shape(cfg["prop"])[0] < np.max(elemprop) + 1
        or np.all(elemprop <= 0)
    ):
        raise ValueError("cfg['prop'] field is missing or insufficient")
    if not isinstance(cfg["e0"], str):
        cfg["e0"] = int(cfg["e0"])
    return cfg


def _unpackdetp(detp, cfg):
    """split the detected photon records into named fields, as mmclab.m"""
    medianum = np.shape(cfg["prop"])[0] - 1
    newdetp = {}
    newdetp["detid"] = detp[0, :].astype(np.int32)
    newdetp["nscat"] = (
        detp[1 : medianum + 1, :].astype(np.int32).T
    )  # number of scattering events per medium
    newdetp["ppath"] = detp[
        medianum + 1 : 2 * medianum + 1, :
    ].T  # partial path lengths per medium
    if cfg.get("ismomentum"):
        newdetp["mom"] = detp[
            2 * medianum + 1 : 3 * medianum + 1, :
        ].T  # momentum transfer per medium
    if cfg.get("issaveexit"):
        newdetp["p"] = detp[-7:-4, :].T  # exit positions
        newdetp["v"] = detp[-4:-1, :].T  # exit directions
    newdetp["w0"] = detp[-1, :].copy()  # the initial packet weight
    if np.ndim(cfg.get("srcpattern")) == (
        4 if cfg.get("srctype") == "pattern3d" else 3
    ):  # several patterns
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
    """Run an MMC simulation with the pre- and post-processing of ``mmclab.m``.

    ``pmmc.run`` only runs the simulation on a prepared mesh. Before running, ``mmclab``:

    * reorients the elements to positive volumes (unless ``cfg['isreoriented']`` is 1), and
      computes ``cfg['facenb']`` and ``cfg['evol']`` when missing;
    * takes ``elemprop`` from the 5th column of ``elem`` when present;
    * finds the initial element ``cfg['e0']``, ray-tracing to the mesh surface for a pencil
      beam outside the mesh, or when ``e0`` is ``'>'``, ``'<'`` or ``'-'``
      (see :func:`mmcraytrace`);
    * for wide-field sources (planar, pattern, fourier, fourierx, fourierx2d; except with
      OptiX) and for sources outside the mesh, adds a source domain to the mesh
      (:func:`mmcsrcdomain`, :func:`mmcaddsrc`).

    ``mmclab.m`` chooses the saved outputs from its number of outputs; here they are chosen
    by ``cfg['issavedet']`` (detected photons), ``cfg['issaveseed']`` (their seeds) and
    ``cfg['debuglevel']`` containing ``'M'`` (trajectories).

    Parameters
    ----------
    cfg : dict, list of dict, or str
        The simulation settings, with the keys of ``pmmc.run`` (see ``help(pmmc.run)``):
        at least ``nphoton``, ``node``, ``elem``, ``elemprop`` (or a 5-column ``elem``),
        ``prop``, ``srcpos``, ``srcdir``, ``tstart``, ``tend`` and ``tstep``. A list runs
        each one. The strings ``'gpuinfo'`` and ``'version'`` return ``pmmc.gpuinfo()`` and
        ``pmmc.version()``.
    option : str, optional
        ``'prep'`` returns the prepared ``cfg`` without running it; ``'sse'`` (or ``'omp'``),
        ``'opencl'``, ``'cuda'`` or ``'optix'`` runs on that backend, the same as setting
        ``cfg['compute']``.

    Returns
    -------
    res : dict or list of dict
        The outputs of ``pmmc.run``, post-processed:

        * ``flux``: the fluence rate (or the quantity of ``cfg['outputtype']``) on the nodes,
          elements or grid; complex for RF simulations (``cfg['omega'] > 0``)
        * ``dref``: the diffuse reflectance on the surface, when ``cfg['issaveref']`` is 1
        * ``jmua``, ``jd``, ...: adjoint Jacobians, ``(nn, nch)`` on the mesh or
          ``(nx, ny, nz, ngate, nch)`` in grid mode, where the ``nch`` channels are all
          source-detector pairs (detector index fastest) or the ``[srcid, detid]`` rows of
          ``cfg['sdmap']``; complex for RF
        * ``detp``: the detected photons, a dict with ``detid``, ``nscat``, ``ppath``,
          ``mom``, ``p``, ``v``, ``w0``, ``prop``, ``unitinmm`` and ``data`` (the raw
          records)
        * ``seeds``: the RNG seeds of the detected photons, one column per photon
        * ``traj``: the trajectories sorted by photon, a dict with ``id``, ``pos``, ``eid``
          and ``data``
        * ``cfg``: the prepared ``cfg``; its ``node``/``elem`` include any added source
          domain, and nodal outputs refer to this mesh

    Raises
    ------
    ValueError
        If a required key is missing, the mesh has degenerate elements, or a wide-field
        source sits on or inside the mesh surface so that no source domain can be added.
    ImportError
        If ``iso2mesh`` is not installed.

    Examples
    --------
    >>> import numpy as np, iso2mesh as i2m, pmmc
    >>> node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 60], 10, 100)
    >>> cfg = {"nphoton": 1000000, "node": node, "elem": elem,
    ...        "elemprop": np.ones(elem.shape[0]), "tstart": 0, "tend": 5e-9, "tstep": 5e-9,
    ...        "srcpos": [30, 30, 0], "srcdir": [0, 0, 1],
    ...        "prop": [[0, 0, 1, 1], [0.005, 1, 0.01, 1.37]]}
    >>> res = pmmc.mmclab(cfg)
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

    if option not in (None, "prep", "omp", "sse", "opencl", "cuda", "optix"):
        raise ValueError("option is not recognized: " + str(option))
    if option not in (None, "prep"):
        # option='sse'/'opencl'/'cuda'/'optix' picks the backend, as cfg['compute'] (mmclab.m);
        # set before _prepare, which skips the source domain for optix
        cfg = dict(cfg, compute=option.replace("omp", "sse"))

    cfg = _prepare(cfg)
    if option == "prep":
        return cfg

    if isinstance(cfg.get("gpuid"), (int, np.integer)) and cfg["gpuid"] < -1:
        cfg["gpuid"] = 1
    res = _run(cfg)

    if (
        "fluximag" in res
    ):  # RF output: return the fluence as a complex array, as mmclab.m does
        res["flux"] = res["flux"] + 1j * np.asarray(
            res.pop("fluximag"), dtype=np.float64
        )
    for key in [
        k for k in res if k.endswith("_im") and k[:-3] in res
    ]:  # RF adjoint Jacobians (jmua_im, jd_im, ...)
        res[key[:-3]] = res[key[:-3]] + 1j * np.asarray(res.pop(key))
    if "detp" in res and cfg.get("issaveexit", 0) != 2 and np.size(res["detp"]) > 0:
        res["detp"] = _unpackdetp(np.asarray(res["detp"]), cfg)
    if "traj" in res and np.size(res["traj"]) > 0:
        res["traj"] = _unpacktraj(res["traj"])
    res["cfg"] = cfg
    return res
