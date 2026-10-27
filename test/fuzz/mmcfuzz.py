#!/usr/bin/env python3
"""
mmcfuzz.py - mutation fuzzer for mmc: mutated JSON inputs and command-line options

usage:
    mmcfuzz.py --bin MMC [options]

Seeds are the test/builtin inputs and the built-in benchmarks (mmc --bench NAME --dumpjson), including
the implicit-MMC ones. Each case mutates a seed field by field (type swaps, extreme numbers, wrong
array shapes, deleted keys) and in mmc-specific ways (MeshROI shapes and values, MeshROIPrim rows,
the ROI medium, sources, output types), adds random command-line flags, and runs it. A case is
reported when mmc crashes, hangs or prints a sanitizer error; an input mmc rejects with an
"MMC ERROR" is a handled case. Reported cases are saved under OUT/caseNNNN with the input, the
arguments and the log.

Build mmc with the sanitizers to catch memory errors, e.g. in src/
    make CC='gcc -fsanitize=address,undefined -fno-omit-frame-pointer'
and fuzz the CPU path (the default -G -1), or a GPU backend with --extra '-c opencl -G 1'.
"""
import argparse, os, json, random, subprocess, glob, shutil, re, math

ap = argparse.ArgumentParser(description="mutation fuzzer for mmc")
ap.add_argument("--bin", required=True, help="mmc binary to fuzz")
ap.add_argument(
    "--builtin",
    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "builtin"),
    help="folder of seed inputs [../builtin]",
)
ap.add_argument(
    "--bench",
    default="dmmc-cube60,dmmc-cube60b,edgeimmc,nodeimmc,faceimmc",
    help="built-in benchmarks added as seeds [dmmc-cube60,dmmc-cube60b,edgeimmc,nodeimmc,faceimmc]",
)
ap.add_argument(
    "--extra", default="-G -1", help="backend options for every run ['-G -1']"
)
ap.add_argument("-n", "--cases", type=int, default=200, help="number of cases [200]")
ap.add_argument("--seed", type=int, default=1, help="fuzzer RNG seed [1]")
ap.add_argument("--out", default="fuzzout", help="output folder [fuzzout]")
ap.add_argument("--timeout", type=float, default=60, help="seconds per run [60]")
opt = ap.parse_args()

rng = random.Random(opt.seed)
os.makedirs(opt.out, exist_ok=True)
EXE = os.path.abspath(opt.bin)
env = dict(
    os.environ,
    ASAN_OPTIONS=os.environ.get("ASAN_OPTIONS", "protect_shadow_gap=0:detect_leaks=0"),
    UBSAN_OPTIONS=os.environ.get("UBSAN_OPTIONS", "print_stacktrace=1"),
    OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"),
)

seeds = [
    json.load(open(f))
    for f in sorted(glob.glob(os.path.join(opt.builtin, "*.json")))
    if os.path.getsize(f) < 2e6 and "atlas" not in f and "colin27" not in f
]

for b in filter(None, opt.bench.split(",")):
    p = subprocess.run(
        [EXE, "--bench", b, "--dumpjson", "-"], capture_output=True, env=env
    )
    try:
        seeds.append(
            json.loads(re.sub(r"\x1b\[[0-9;]*m", "", p.stdout.decode("latin1")))
        )
    except ValueError:
        print("can not load benchmark", b)

ODD = [
    0,
    -1,
    1,
    1e-30,
    -1e-30,
    1e30,
    -1e30,
    2**31,
    2**32 + 1,
    2**63,
    0.5,
    1e-6,
    1e6,
    3,
    float("nan"),
    float("inf"),
]
OTYPES = [
    "flux",
    "fluence",
    "energy",
    "jacobian",
    "wl",
    "wp",
    "rf",
    "rfmus",
    "x",
    "f",
    "e",
    "garbage",
    7,
]
SRCTYPES = [
    "pencil",
    "isotropic",
    "cone",
    "gaussian",
    "planar",
    "pattern",
    "fourier",
    "arcsine",
    "disk",
    "fourierx",
    "fourierx2d",
    "zgaussian",
    "line",
    "slit",
    "pencilarray",
    "pattern3d",
    "hyperboloid",
    "ring",
    "nosuch",
]
FLAGS = [
    ["-M", ["G", "S", "P", "H", "B", "X"]],
    ["-C", ["0", "1", "2"]],
    ["-b", ["0", "1", "2", "3"]],
    ["-O", list("FELXJPWRQAVDZ")],
    ["-D", ["T", "M", "TPM", "A", "X"]],
    ["-x", ["0", "1"]],
    ["-X", ["0", "1"]],
    ["-d", ["0", "1"]],
    ["-e", ["0", "1e-9", "0.5"]],
    ["-U", ["0", "1"]],
    ["-q", ["0", "1"]],
    ["-V", ["0", "1", "2"]],
    ["-k", ["0", "1"]],
    ["-m", ["0", "1"]],
    ["-u", ["1", "0.005", "-1"]],
    ["-S", ["0", "1"]],
]


def paths(o, p=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from paths(v, p + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o[:8]):
            yield from paths(v, p + (i,))
    yield p


def get(o, p):
    for k in p:
        o = o[k]
    return o


def meshobj(d):
    """the mesh object of an input: Mesh, or Shapes as written by --dumpjson"""
    for k in ("Mesh", "Shapes"):
        if isinstance(d.get(k), dict):
            return d[k]
    return None


def mutate_roi(d):
    """implicit MMC: MeshROI shapes and values, MeshROIPrim rows and the ROI medium"""
    m = meshobj(d)

    if m is None:
        return

    r = rng.random()

    if r < 0.3:
        rows, cols = rng.choice([1, 2, 6, 12, 16, 1000]), rng.choice([1, 4, 6, 3, 7])
        m["MeshROI"] = [
            [
                rng.choice([0, 0, 0.1, 0.5, -1, 1e30, 0.0, float("nan")])
                for _ in range(cols)
            ]
            for _ in range(rows)
        ]
    elif r < 0.6:
        rows = []

        for _ in range(rng.randint(1, 5)):
            a = [rng.uniform(-1, 2) for _ in range(3)]
            b = a[:] if rng.random() < 0.3 else [rng.uniform(-1, 2) for _ in range(3)]
            rows.append(
                a + b + [rng.choice([0.05, 0.3, 0, -0.1, 1e6, 1e-9, float("nan")])]
            )

        if rng.random() < 0.2:
            rows.append(rng.choice([[1, 2, 3], [], "x", [1] * 8]))

        m["MeshROIPrim"] = rows
    elif r < 0.75:
        m.pop("MeshROI", None)

    dom = d.get("Domain") if isinstance(d.get("Domain"), dict) else None

    if (
        dom
        and isinstance(dom.get("Media"), list)
        and dom["Media"]
        and isinstance(dom["Media"][-1], dict)
        and rng.random() < 0.5
    ):
        dom["Media"][-1]["n"] = rng.choice([1.0, 1.37, 1.5, 3.0, 0.5, 0, -1])


def mutate_json(d):
    for _ in range(rng.randint(1, 3)):
        ps = [p for p in paths(d) if p]

        if not ps:
            break

        p = rng.choice(ps)

        try:
            v, par, r = get(d, p), get(d, p[:-1]), rng.random()

            if r < 0.3 and isinstance(v, (int, float)) and not isinstance(v, bool):
                par[p[-1]] = (
                    rng.choice(ODD)
                    if rng.random() < 0.5
                    else v * rng.choice([0, -1, 2, 10, 0.1, 1e3])
                )
            elif r < 0.45:
                par[p[-1]] = rng.choice(["abc", "", [], {}, None, True, [1, 2], -5])
            elif r < 0.6 and isinstance(v, list):
                par[p[-1]] = (
                    v[: rng.randint(0, len(v))]
                    if v and rng.random() < 0.5
                    else v + [rng.choice(ODD)]
                )
            elif r < 0.7 and isinstance(par, dict):
                del par[p[-1]]
        except Exception:
            pass

    if rng.random() < 0.4:
        mutate_roi(d)

    sess = d.get("Session") if isinstance(d.get("Session"), dict) else None
    optode = d.get("Optode") if isinstance(d.get("Optode"), dict) else {}
    src = optode.get("Source") if isinstance(optode.get("Source"), dict) else None

    if sess is not None and rng.random() < 0.3:
        sess["OutputType"] = rng.choice(OTYPES)

    if src is not None and rng.random() < 0.2:
        src["Type"] = rng.choice(SRCTYPES)

    if isinstance(d.get("Forward"), dict) and rng.random() < 0.15:
        d["Forward"].update(
            rng.choice(
                [{"T1": 1e-6}, {"Dt": 1e-12}, {"T1": 5e-11, "Dt": 1e-11}, {"T0": 1e-9}]
            )
        )

    return d


def mutate_cli():
    a = []

    for _ in range(rng.randint(0, 4)):
        f, vals = rng.choice(FLAGS)
        a += [f, rng.choice(vals)]

    return a


def run(args, wd):
    try:
        p = subprocess.run(
            [EXE] + args, cwd=wd, env=env, capture_output=True, timeout=opt.timeout
        )
        txt, rc = (p.stdout + p.stderr).decode("latin1"), p.returncode
    except subprocess.TimeoutExpired:
        txt, rc = "TIMEOUT", "timeout"

    txt = re.sub(r"\x1b\[[0-9;]*m", "", txt)
    return {
        "rc": rc,
        "san": [
            l for l in txt.splitlines() if "Sanitizer" in l or "runtime error" in l
        ],
        "err": re.findall(r"MMC ERROR\(?(-?\d+)\)?", txt),
        "log": txt,
    }


stats = {}

for c in range(opt.cases):
    si = rng.randrange(len(seeds))
    d = mutate_json(json.loads(json.dumps(seeds[si])))
    args = (
        [
            "-f",
            "in.json",
            "-n",
            rng.choice(["1e2", "1e3", "1e4"]),
            "-E",
            "1648335518",
            "-F",
            "bin",
            "-s",
            "o",
        ]
        + opt.extra.split()
        + mutate_cli()
    )
    wd = os.path.join(opt.out, "w")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)

    with open(os.path.join(wd, "in.json"), "w") as fp:
        json.dump(d, fp, allow_nan=True)

    r = run(args, wd)
    found = []

    if r["rc"] == "timeout":
        found.append("hang")
    elif isinstance(r["rc"], int) and r["rc"] < 0:
        found.append("crash: signal %d" % -r["rc"])

    if r["san"]:
        found.append("sanitizer: " + r["san"][0][:120])

    key = found[0].split(":")[0] if found else ("rejected" if r["err"] else "ok")
    stats[key] = stats.get(key, 0) + 1

    if found:
        cd = os.path.join(opt.out, "case%04d" % c)
        os.makedirs(cd, exist_ok=True)
        shutil.copy(os.path.join(wd, "in.json"), cd)
        open(os.path.join(cd, "args.txt"), "w").write(
            " ".join(args) + "\n" + "\n".join(found) + "\n"
        )
        open(os.path.join(cd, "mmc.log"), "w").write(r["log"][-6000:])

    print("%4d seed%-3d %-8s %s" % (c, si, key, "; ".join(found)), flush=True)

print("SUMMARY", json.dumps(stats))
