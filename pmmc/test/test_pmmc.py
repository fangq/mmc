"""Minimal pmmc tests on the CPU (SSE/OpenMP) path, no GPU needed; the GPU
backends (opencl, and cuda when built with CMAKE_ARGS="-DBUILD_CUDA=ON") are
tested on device 1 when usable, and skipped otherwise

Run with: python3 -m unittest discover -s pmmc/test
"""

import os
import unittest
from itertools import permutations

os.environ.setdefault("OMP_NUM_THREADS", "1")  # single thread: reproducible

import numpy as np  # noqa: E402

try:
    import pmmc  # noqa: E402
except ImportError:
    pmmc = None


def box_mesh(size=20.0):
    """a cube [0,size]^3 split into 6 Kuhn tetrahedra"""
    node = (
        np.array(
            [[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)],
            dtype=np.float32,
        )
        * size
    )
    elem = []

    for perm in permutations(range(3)):
        v = [0, 0, 0]
        tet = [0]

        for ax in perm:
            v[ax] = 1
            tet.append(v[0] * 4 + v[1] * 2 + v[2])

        elem.append(tet)

    return node, np.asarray(elem, dtype=np.int32) + 1


def base_cfg():
    node, elem = box_mesh()
    return {
        "nphoton": 20000,
        "seed": 12345,
        "node": node,
        "elem": elem,
        "elemprop": np.ones(elem.shape[0], dtype=np.int32),
        "srcpos": np.array([10, 10, 0], dtype=np.float32),
        "srcdir": np.array([0, 0, 1], dtype=np.float32),
        "tstart": 0,
        "tstep": 5e-9,
        "tend": 5e-9,
        "prop": np.array([[0, 0, 1, 1], [0.05, 2, 0.8, 1]], dtype=np.float32),
        "e0": 1,
        "compute": "sse",
    }


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestPMMC(unittest.TestCase):
    def test_flux_is_reproducible(self):
        cfg = base_cfg()
        cfg["nphoton"] = 2000
        a = pmmc.run(cfg)
        b = pmmc.run(cfg)
        self.assertIn("flux", a)
        self.assertTrue(np.all(a["flux"] >= 0))
        self.assertTrue(np.array_equal(a["flux"], b["flux"]))

    def test_absorbed_fraction(self):
        # 1 - mean exit weight of the same box simulated by mcx is 0.574 (+/-0.005)
        cfg = base_cfg()
        cfg["outputtype"] = "energy"
        absorbed = pmmc.run(cfg)["flux"].sum()
        self.assertGreater(absorbed, 0.54)
        self.assertLess(absorbed, 0.60)


@unittest.skipIf(pmmc is None, "pmmc is not installed")
class TestMMCLAB(unittest.TestCase):
    def test_option_sets_compute(self):
        """mmclab(cfg, 'sse'/'omp'/'opencl'/...) runs on that backend, as cfg['compute'] (mmclab.m)"""
        try:
            import iso2mesh  # noqa: F401
        except ImportError:
            self.skipTest("pmmc.mmclab needs iso2mesh")

        import sys
        from unittest.mock import patch

        mod = sys.modules["pmmc.mmclab"]  # pmmc.mmclab is the function
        seen = []

        def fakerun(cfg):
            seen.append(cfg["compute"])
            return {"flux": np.zeros(1)}

        cfg = base_cfg()
        cfg.update({"compute": "opencl", "nphoton": 100})
        with patch.object(mod, "_run", fakerun):
            for option in ("sse", "omp", "cuda"):
                pmmc.mmclab(cfg, option)
            pmmc.mmclab(cfg)

        self.assertEqual(seen, ["sse", "sse", "cuda", "opencl"])
        self.assertEqual(cfg["compute"], "opencl")  # the caller's cfg is unchanged
        with self.assertRaises(ValueError):
            pmmc.mmclab(cfg, "gpu")


def has_backend(name):
    """True if pmmc can run on the given GPU backend on device 1"""
    if pmmc is None:
        return False

    cfg = base_cfg()
    cfg.update({"nphoton": 100, "compute": name, "gpuid": 1})

    try:
        return "flux" in pmmc.run(cfg)
    except Exception:
        return False


class TestPMMCGPU(unittest.TestCase):
    def check_absorbed(self, backend):
        if not has_backend(backend):
            self.skipTest("pmmc %s backend is not available" % backend)

        cfg = base_cfg()
        cfg.update(
            {"outputtype": "energy", "compute": backend, "gpuid": 1, "nphoton": 100000}
        )
        absorbed = pmmc.run(cfg)["flux"].sum()
        self.assertGreater(absorbed, 0.56)
        self.assertLess(absorbed, 0.59)

    def test_absorbed_fraction_opencl(self):
        self.check_absorbed("opencl")

    def test_absorbed_fraction_cuda(self):
        self.check_absorbed("cuda")

    def test_sdmap_opencl(self):
        # cfg['sdmap'] rows [srcid, detid] (1-based) select the adjoint Jacobian channels
        if not has_backend("opencl"):
            self.skipTest("pmmc opencl backend is not available")

        cfg = base_cfg()
        cfg.update(
            {
                "compute": "opencl",
                "gpuid": 1,
                "nphoton": 200000,
                "outputtype": "adjoint",
                "basisorder": 1,
                "detpos": np.array(
                    [[5, 10, 0, 2], [15, 10, 0, 2], [10, 15, 0, 2]], dtype=np.float32
                ),
                "detdir": np.tile(np.array([0, 0, 1, 0], dtype=np.float32), (3, 1)),
            }
        )
        full = pmmc.run(cfg)["jmua"]
        cfg["sdmap"] = np.array([[1, 3], [1, 1]])
        sub = pmmc.run(cfg)["jmua"]
        self.assertEqual(full.shape[-1], 3)
        self.assertEqual(sub.shape[-1], 2)

        for c, k in enumerate([2, 0]):
            # the two runs differ only by GPU run-to-run noise
            err = np.linalg.norm(full[..., k] - sub[..., c]) / np.linalg.norm(
                full[..., k]
            )
            self.assertLess(err, 0.05)

        cfg["sdmap"] = np.array([[2, 1]])

        with self.assertRaises(Exception):
            pmmc.run(cfg)


if __name__ == "__main__":
    unittest.main()
