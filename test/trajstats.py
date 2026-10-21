#!/usr/bin/env python3
"""
trajstats.py - statistically compare photon trajectories of two solvers

Usage:
    trajstats.py TRAJ_A TRAJ_B [--sigma 4]
    trajstats.py TRAJ_A --save-stats STATS.json [--solver mcx]

TRAJ_A/B are trajectory files (session_traj.jdt from -D M), solver-neutral
golden files readable by trajtest.py, or statistics files written by
--save-stats (count, mean and variance of each quantity), for example MMC on a box mesh and MCX on
the matching voxel grid. MMC and MCX consume random numbers differently, so
the individual paths differ; instead, the per-photon total path length, the
remaining weight at the end of the path and the exit position are compared as
distributions: each mean must agree within SIGMA combined standard errors.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trajtest import load_traj  # noqa: E402


def stats(photons):
    plen = np.array([np.linalg.norm(np.diff(ph[:, :3], axis=0), axis=1).sum() for ph in photons])
    wend = np.array([ph[-1, 3] for ph in photons])
    pend = np.array([ph[-1, :3] for ph in photons])
    return {"path length": plen, "final weight": wend, "exit x": pend[:, 0], "exit y": pend[:, 1], "exit z": pend[:, 2]}


def summarize(fname):
    """count/mean/var of each quantity, from a trajectory or a saved stats file"""
    if fname.endswith(".json"):
        with open(fname) as fp:
            data = json.load(fp)

        if "stats" in data:
            return data["stats"]

    return {k: {"n": len(v), "mean": float(v.mean()), "var": float(v.var())} for k, v in stats(load_traj(fname)).items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b", nargs="?")
    ap.add_argument("--sigma", type=float, default=4.0)
    ap.add_argument("--save-stats", dest="savestats")
    ap.add_argument("--solver", default="mcx")
    opt = ap.parse_args()

    sa = summarize(opt.a)

    if opt.savestats:
        with open(opt.savestats, "w") as fp:
            json.dump({"solver": opt.solver, "stats": sa}, fp, indent=1)
        return

    sb = summarize(opt.b)
    fail = 0

    for key in sa:
        a, b = sa[key], sb[key]
        se = np.sqrt(a["var"] / a["n"] + b["var"] / b["n"])
        z = abs(a["mean"] - b["mean"]) / max(se, 1e-12)
        ok = z <= opt.sigma
        fail += not ok
        print("  %-13s %10.4f vs %10.4f  (%.2f sigma) %s" % (key, a["mean"], b["mean"], z, "ok" if ok else "FAIL"))

    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
