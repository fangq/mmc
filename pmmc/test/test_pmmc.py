"""Minimal pmmc tests on the CPU (SSE/OpenMP) path, no GPU needed

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
    node = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=np.float32) * size
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


if __name__ == "__main__":
    unittest.main()
