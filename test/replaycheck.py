#!/usr/bin/env python3
"""
replaycheck.py - check the weights of the replayed photons of the mmc binary

Usage:
    replaycheck.py MMC WORKDIR [extra mmc parameters, e.g. -G -1 or -c opencl -G 1]

A two-layer box (traj_cases/twolayer.json tessellated by mkmeshcase.py, mua 0.02/mm
over 0.05/mm) with a detector is simulated once to save the detected photons and their
seeds (.mch), then replayed with -O L (wl) and -O J (jacobian) without normalization.
Every replayed photon i deposits w_i*L along its path, so the total of the map must
equal sum_i w_i*sum_m L_im, where L_im are its partial paths (mesh units) and
w_i = w0_i*exp(-sum_m mua_m*L_im*unitinmm) its detected weight, both read from the
.mch file. jacobian (the weighted path length, as in mcx) must equal wl.

This catches the replay reading the scattering-count columns of the .mch file as the
partial paths (it needs 2 or more media), and applying unitinmm twice to the weights;
the case is run at LengthUnit 1 and 2 (the mesh in 2 mm units).

A time-resolved wl replay with a time window starting after 0 (4 gates) must put
every photon in the gate of its time of flight t_i = sum_m n_m*L_im*unitinmm/c minus
the start of the window, and drop the photons detected outside of the window.

On the CPU (-G -1), this also catches a replayed photon not following the path of the
detected one, e.g. when it uses random numbers buffered by the previous photon.

The exit code is the number of failed checks.
"""

import json
import os
import subprocess
import sys

import numpy as np

HEADER = np.dtype(
    [
        ("magic", "S4"), ("version", "<u4"), ("maxmedia", "<u4"), ("detnum", "<u4"),
        ("colcount", "<u4"), ("totalphoton", "<u4"), ("detected", "<u4"),
        ("savedphoton", "<u4"), ("unitinmm", "<f4"), ("seedbyte", "<u4"),
        ("normalizer", "<f4"), ("srcnum", "<u4"), ("respin", "<i4"),
        ("savedetflag", "<u4"), ("reserved", "<i4", 2),
    ]
)
RTOL = 1e-3  # float32 accumulation over a few hundred photons
C0 = 299792458e3  # speed of light in vacuum, mm/s


def loadmch(fname):
    """the header and the detected photon records (one row per photon) of an .mch file"""
    with open(fname, "rb") as f:
        his = np.frombuffer(f.read(HEADER.itemsize), HEADER)[0]
        data = np.frombuffer(f.read(4 * his["colcount"] * his["savedphoton"]), "<f4")
    return his, data.reshape(his["savedphoton"], his["colcount"])


def run(mmc, workdir, args):
    res = subprocess.run([mmc] + args, cwd=workdir, capture_output=True, text=True)
    if res.returncode:
        sys.exit("mmc %s failed:\n%s%s" % (" ".join(args), res.stdout, res.stderr))


def main():
    mmc, workdir, extra = sys.argv[1], sys.argv[2], sys.argv[3:]
    testdir = os.path.dirname(os.path.abspath(__file__))
    tag = "replay" + "".join(a.strip("-") for a in extra)
    meshcase = os.path.join(workdir, tag + "_mesh.json")
    subprocess.run(
        [sys.executable, os.path.join(testdir, "mkmeshcase.py"),
         os.path.join(testdir, "traj_cases", "twolayer.json"), meshcase],
        check=True,
    )

    with open(meshcase) as f:
        base = json.load(f)

    fail = 0

    for unit in (1, 2):
        cfg = json.loads(json.dumps(base))
        cfg["Session"].update({"ID": "%s_u%d" % (tag, unit), "Photons": 200000})
        cfg["Domain"]["LengthUnit"] = unit
        cfg["Mesh"]["MeshNode"] = (np.asarray(cfg["Mesh"]["MeshNode"], float) / unit).tolist()
        cfg["Optode"]["Source"]["Pos"] = [10.0 / unit, 10.0 / unit, 0]
        cfg["Optode"]["Detector"] = [{"Pos": [10.0 / unit, 15.0 / unit, 0], "R": 2.0 / unit}]
        inp = "%s_u%d.json" % (tag, unit)

        with open(os.path.join(workdir, inp), "w") as f:
            json.dump(cfg, f)

        sess = cfg["Session"]["ID"]
        run(mmc, workdir, ["-f", inp, "-d", "1", "-q", "1", "-S", "0", "-s", sess + "_base"] + extra)
        his, det = loadmch(os.path.join(workdir, sess + "_base.mch"))
        nmed = int(his["maxmedia"])

        if nmed != 2 or his["savedphoton"] < 50:
            print("unexpected baseline: %d media, %d detected photons" % (nmed, his["savedphoton"]))
            fail += 1
            continue

        # columns: detid, nscat x nmed, ppath x nmed, ..., w0
        ppath = det[:, 1 + nmed:1 + 2 * nmed].astype(float)
        mua = np.array([m["mua"] for m in cfg["Domain"]["Media"][1:]])
        nidx = np.array([m["n"] for m in cfg["Domain"]["Media"][1:]])
        weight = det[:, -1] * np.exp(-(ppath * unit) @ mua)
        expected = (weight[:, None] * ppath).sum()
        tof = (ppath * unit) @ nidx / C0

        total = {}

        for otype in ("L", "J"):
            out = "%s_%s" % (sess, otype)
            run(mmc, workdir, ["-f", inp, "-E", sess + "_base.mch", "-O", otype, "-U", "0",
                               "-F", "bin", "-s", out] + extra)
            total[otype] = np.fromfile(os.path.join(workdir, out + ".bin")).sum()

        for otype, name in (("L", "wl"), ("J", "jacobian")):
            ok = abs(total[otype] / expected - 1) < RTOL
            print("unitinmm %d, %s: total of the replay map %.6g, sum of w*L %.6g (%d photons) ... %s"
                  % (unit, name, total[otype], expected, len(weight), "ok" if ok else "FAIL"))
            fail += not ok

        # time-resolved: the window starts at the 30th percentile of the times of flight
        ngate = 4
        t0 = float(np.quantile(tof, 0.3))
        dt = (float(np.quantile(tof, 0.95)) - t0) / ngate
        gate = np.minimum(np.floor((tof - t0) / dt), ngate - 1)
        inwin = (tof >= t0) & (tof <= t0 + ngate * dt)
        wl = (weight[:, None] * ppath).sum(1)
        gexpected = np.array([wl[inwin & (gate == g)].sum() for g in range(ngate)])
        cfg["Forward"] = {"T0": t0, "T1": t0 + ngate * dt, "Dt": dt}
        tinp = "%s_u%d_t.json" % (tag, unit)

        with open(os.path.join(workdir, tinp), "w") as f:
            json.dump(cfg, f)

        out = "%s_Lt" % sess
        run(mmc, workdir, ["-f", tinp, "-E", sess + "_base.mch", "-O", "L", "-U", "0",
                           "-F", "bin", "-s", out] + extra)
        gtotal = np.fromfile(os.path.join(workdir, out + ".bin")).reshape(ngate, -1).sum(1)
        ok = np.allclose(gtotal, gexpected, rtol=RTOL, atol=RTOL * gexpected.max())
        print("unitinmm %d, wl in %d gates from %.3g ns: %s, expected %s ... %s"
              % (unit, ngate, t0 * 1e9, np.array2string(gtotal, precision=5),
                 np.array2string(gexpected, precision=5), "ok" if ok else "FAIL"))
        fail += not ok

    return fail


if __name__ == "__main__":
    sys.exit(main())
