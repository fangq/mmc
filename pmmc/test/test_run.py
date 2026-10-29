"""Unit tests of the pmmc binding interface (run, gpuinfo, version and the mmclab
front end): the cfg parsing, the errors, the output keys, shapes and types, and short
absorbed-fraction end-to-end checks. Everything runs on the CPU (compute='sse') with a
few thousand photons, so no GPU is needed; the mmclab tests need iso2mesh.

Run with: python3 -m unittest discover -s pmmc/test
"""

import os
import subprocess
import sys
import unittest

os.environ.setdefault("OMP_NUM_THREADS", "1")  # single thread: reproducible

import numpy as np  # noqa: E402

try:
    import pmmc  # noqa: E402
except ImportError:
    pmmc = None

try:
    import iso2mesh  # noqa: E402
except ImportError:
    iso2mesh = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_pmmc import base_cfg  # noqa: E402  the 20 mm cube of 6 tetrahedra

PROP = np.array([[0, 0, 1, 1], [0.05, 2, 0.8, 1]], np.float32)  # as base_cfg
NULLDEV = "NUL" if sys.platform.startswith("win") else "/dev/null"


def cfg(**kw):
    """base_cfg() with a few thousand photons and a silent log; keywords override"""
    c = base_cfg()
    c.update(nphoton=2000, flog=NULLDEV)
    c.update(kw)
    return c


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestInfo(unittest.TestCase):
    def test_version(self):
        ver = pmmc.version()
        self.assertIsInstance(ver, str)
        self.assertTrue(ver.startswith("v"), ver)

    def test_gpuinfo_is_a_list(self):
        try:
            info = pmmc.gpuinfo()
        except Exception:  # no OpenCL runtime
            self.skipTest("no OpenCL runtime")
        self.assertIsInstance(info, list)
        if info:
            self.assertIn("name", info[0])

    def test_docstrings(self):
        for name in ("run", "gpuinfo", "version"):
            self.assertIn("Returns", getattr(pmmc, name).__doc__ or "", name)


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestErrors(unittest.TestCase):
    """invalid settings raise Python exceptions, not crashes or silent runs"""

    def test_missing_mesh(self):
        # a missing node or elem used to crash the interpreter
        for key in ("node", "elem", "elemprop"):
            c = cfg()
            del c[key]
            with self.subTest(missing=key):
                with self.assertRaises(ValueError):
                    pmmc.run(c)

    def test_missing_prop(self):
        c = cfg()
        del c["prop"]
        with self.assertRaises((ValueError, RuntimeError)):
            pmmc.run(c)

    def test_bad_values(self):
        for kw in (
            dict(method="xyz"),
            dict(compute="xyz"),
            dict(debuglevel="Q"),
            dict(outputtype="xyz"),
            dict(srctype="xyz"),
        ):
            with self.subTest(kw=list(kw)):
                with self.assertRaises(ValueError):
                    pmmc.run(cfg(**kw))


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestInterface(unittest.TestCase):
    def test_dict_and_keywords(self):
        """run(cfg) and run(**cfg) are the same simulation"""
        np.testing.assert_array_equal(
            pmmc.run(cfg())["flux"], pmmc.run(**cfg())["flux"]
        )

    def test_cfg_unchanged(self):
        c = cfg(detpos=[[10, 8, 0, 2]])
        before = {
            k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v)
            for k, v in c.items()
        }
        pmmc.run(c)
        self.assertEqual(set(c), set(before))
        for k in c:
            np.testing.assert_array_equal(c[k], before[k], err_msg=k)

    def test_computed_mesh_data(self):
        """facenb, evol and e0 are computed when absent"""
        c = cfg()
        del c["e0"]
        np.testing.assert_allclose(pmmc.run(c)["flux"], pmmc.run(cfg())["flux"])

    def test_output_shapes(self):
        """'grid' gives voxels covering the mesh; the mesh tracers give nodal or elemental values"""
        self.assertEqual(pmmc.run(cfg())["flux"].shape, (21, 21, 21, 1))
        self.assertEqual(pmmc.run(cfg(tstep=1e-9))["flux"].shape, (21, 21, 21, 5))
        self.assertEqual(pmmc.run(cfg(method="elem"))["flux"].shape, (8, 1))
        self.assertEqual(
            pmmc.run(cfg(method="elem", basisorder=0))["flux"].shape, (6, 1)
        )
        self.assertEqual(pmmc.run(cfg())["flux"].dtype, np.float64)

    def test_energy_is_the_same_on_grid_and_mesh(self):
        """the absorbed energy does not depend on where it is accumulated"""
        grid = pmmc.run(cfg(outputtype="energy", nphoton=20000))["flux"].sum()
        elem = pmmc.run(
            cfg(outputtype="energy", nphoton=20000, method="elem", basisorder=0)
        )["flux"].sum()
        self.assertAlmostEqual(grid, elem, delta=0.005)

    def test_detected_photons(self):
        """detpos alone saves the detected photons; the rows follow the enabled fields"""
        res = pmmc.run(
            cfg(detpos=[[10, 8, 0, 2]], issaveseed=1, ismomentum=1, issaveexit=1)
        )
        ndet = res["detp"].shape[1]
        self.assertGreater(ndet, 0)
        self.assertEqual(
            res["detp"].shape[0], 1 + 3 * 1 + 3 + 3 + 1
        )  # detid, nscat/ppath/mom, p, v, w0
        self.assertEqual(res["seeds"].dtype, np.uint8)
        self.assertEqual(res["seeds"].shape[1], ndet)
        np.testing.assert_array_equal(res["detp"][0], 1)
        self.assertEqual(
            pmmc.run(cfg(detpos=[[10, 8, 0, 2]]))["detp"].shape[0], 4
        )  # detid, nscat, ppath, w0
        self.assertNotIn("detp", pmmc.run(cfg(detpos=[[10, 8, 0, 2]], issavedet=0)))

    def test_no_fluence(self):
        res = pmmc.run(cfg(issave2pt=0, detpos=[[10, 8, 0, 2]]))
        self.assertNotIn("flux", res)
        self.assertIn("detp", res)

    def test_trajectories(self):
        res = pmmc.run(cfg(nphoton=100, debuglevel="M"))
        traj = res["traj"]
        self.assertEqual(traj.shape[0], 6)  # id, x, y, z, weight, eid
        pos = traj[1:4]
        self.assertTrue(np.all(pos >= -1e-3) and np.all(pos <= 20 + 1e-3))
        eid = traj[5].view(
            np.int32
        )  # the enclosing element, negative where the photon exits
        self.assertTrue(np.all(((eid >= 1) & (eid <= 6)) | (eid < 0)))
        self.assertTrue(np.any(eid > 0))

    def test_reflectance(self):
        res = pmmc.run(cfg(issaveref=1, method="elem"))
        self.assertEqual(res["dref"].shape, (12, 1))  # 12 surface triangles of the cube
        self.assertGreater(res["dref"].sum(), 0)

    def test_rf_parts(self):
        res = pmmc.run(cfg(omega=2 * np.pi * 1e8, method="elem"))
        self.assertEqual(res["flux"].shape, res["fluximag"].shape)
        self.assertTrue(np.any(res["fluximag"] != 0))

    def test_replay(self):
        base = pmmc.run(cfg(nphoton=20000, detpos=[[10, 8, 0, 2]], issaveseed=1))
        jac = pmmc.run(
            cfg(
                seed=base["seeds"],
                detphotons=base["detp"],
                outputtype="jacobian",
                detpos=[[10, 8, 0, 2]],
            )
        )
        self.assertTrue(np.all(np.isfinite(jac["flux"])))
        self.assertGreater(np.abs(jac["flux"]).sum(), 0)

    def test_flog_null_device_is_silent(self):
        script = (
            "import sys; sys.path[:0] = %r; from test_pmmc import base_cfg; import pmmc; "
            "c = base_cfg(); c.update(nphoton=200, flog=%r); pmmc.run(c)"
            % (sys.path, NULLDEV)
        )
        out = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, timeout=300
        )
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.strip(), "")


@unittest.skipIf(pmmc is None or iso2mesh is None, "pmmc.mmclab needs iso2mesh")
class TestMmclab(unittest.TestCase):
    """the mesh preparation and the post-processing of the mmclab front end"""

    @classmethod
    def setUpClass(cls):
        node, face, elem = iso2mesh.meshabox([0, 0, 0], [20, 20, 20], 10)
        cls.mesh = dict(node=node, elem=elem[:, :4], elemprop=np.ones(elem.shape[0]))

    def mcfg(self, **kw):
        c = dict(
            nphoton=2000,
            srcpos=[10, 10, 0],
            srcdir=[0, 0, 1],
            tstart=0,
            tend=5e-9,
            tstep=5e-9,
            prop=PROP,
            flog=NULLDEV,
            compute="sse",
            seed=1,
        )
        c.update(self.mesh)
        c.update(kw)
        return c

    def test_prep(self):
        prep = pmmc.mmclab(self.mcfg(), "prep")
        ne = self.mesh["elem"].shape[0]
        self.assertEqual(prep["isreoriented"], 1)
        self.assertEqual(np.shape(prep["facenb"]), (ne, 4))
        self.assertEqual(np.shape(prep["evol"]), (ne,))
        self.assertTrue(np.all(np.asarray(prep["evol"]) > 0))
        self.assertTrue(1 <= prep["e0"] <= ne)

    def test_named_outputs(self):
        res = pmmc.mmclab(
            self.mcfg(detpos=[[10, 8, 0, 2]], issaveseed=1, debuglevel="M")
        )
        detp = res["detp"]
        ndet = detp["detid"].shape[0]
        self.assertGreater(ndet, 0)
        self.assertEqual(detp["ppath"].shape, (ndet, 1))
        self.assertEqual(res["seeds"].shape[1], ndet)
        self.assertTrue(np.all(np.diff(res["traj"]["id"].astype(np.int64)) >= 0))
        self.assertEqual(res["traj"]["pos"].shape[1], 3)
        self.assertIn("cfg", res)

    def test_rf_is_complex(self):
        res = pmmc.mmclab(self.mcfg(omega=2 * np.pi * 1e8, method="elem"))
        self.assertTrue(np.iscomplexobj(res["flux"]))
        self.assertNotIn("fluximag", res)

    def test_widefield_source_domain(self):
        """a planar source outside the mesh gets a source domain (label -1) added"""
        src = dict(
            srctype="planar",
            srcpos=[5, 5, -2],
            srcparam1=[10, 0, 0, 0],
            srcparam2=[0, 10, 0, 0],
        )
        res = pmmc.mmclab(self.mcfg(**src))
        elemprop = np.ravel(res["cfg"]["elemprop"])
        self.assertGreater(np.sum(elemprop < 0), 0)
        self.assertGreater(res["cfg"]["elem"].shape[0], self.mesh["elem"].shape[0])
        self.assertGreater(np.nansum(res["flux"]), 0)

    def test_widefield_source_on_surface(self):
        """a planar source on the mesh surface has no room for a source domain: an error"""
        src = dict(
            srctype="planar",
            srcpos=[5, 5, 0],
            srcparam1=[10, 0, 0, 0],
            srcparam2=[0, 10, 0, 0],
        )
        with self.assertRaises(ValueError):
            pmmc.mmclab(self.mcfg(**src))

    def test_info_commands(self):
        self.assertTrue(pmmc.mmclab("version").startswith("v"))


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestEndToEnd(unittest.TestCase):
    """short absorbed-fraction checks (the cube itself is checked in test_pmmc)"""

    def test_boundary_reflection(self):
        """Fresnel reflection keeps more light inside; isreflect=0 is an index-matched surface"""
        p137 = [[0, 0, 1, 1], [0.05, 2, 0.8, 1.37]]
        reflect = pmmc.run(cfg(outputtype="energy", nphoton=20000, prop=p137))[
            "flux"
        ].sum()
        noreflect = pmmc.run(
            cfg(outputtype="energy", nphoton=20000, prop=p137, isreflect=0)
        )["flux"].sum()
        matched = pmmc.run(cfg(outputtype="energy", nphoton=20000))[
            "flux"
        ].sum()  # n=1 inside and out
        self.assertGreater(reflect, noreflect + 0.05)
        self.assertAlmostEqual(noreflect, matched, delta=1e-6)  # the same photon paths
        self.assertGreater(matched, 0.54)  # 0.574 from mcx, see test_pmmc
        self.assertLess(matched, 0.60)


if __name__ == "__main__":
    unittest.main()
