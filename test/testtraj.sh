#!/bin/sh
#
# testtraj.sh - step-by-step photon trajectory tests for MMC
#
# usage: ./testtraj.sh [-u] [extra mmc parameters]
#
# each box-shaped case under traj_cases/ (shared with mcx/test/traj_cases) is
# tessellated into a tetrahedral mesh by mkmeshcase.py and simulated on the
# CPU with a single thread (OMP_NUM_THREADS=1) and a fixed seed, so that the
# saved trajectory (-D M) is deterministic. Two kinds of checks are run:
#
#   1. the 20-photon trajectory is checked against physical invariants and
#      compared step-by-step with the MMC golden file golden/<case>.json;
#   2. a larger run is compared statistically with MCX on the matching voxel
#      grid (golden/<case>_mcxstats.json): mean path length, final weight and
#      exit position must agree within SIGMA standard errors. MMC consumes
#      random numbers differently from MCX, so paths only agree in distribution.
#
# -u regenerates the MMC golden files.
#
# environment variables:
#    MMC:     path to the mmc binary (default ../bin/mmc)
#    MMCOPT:  device options (default "-G -1", the SSE/OpenMP CPU path)
#    PYTHON:  python interpreter with numpy (default python3)
#    REFTOL:  position (mm)/relative weight tolerance against the golden data (default 1e-3)
#    MAXFORK: maximum fraction of photons allowed to fork from the golden data
#             (default 0, i.e. the same backend must reproduce every step; use
#             e.g. REFTOL=5e-2 MAXFORK=0.1 to compare another backend, see trajtest.py)
#    GOLDENTAG: suffix of the MMC golden files, e.g. "_opencl" for the OpenCL
#             path, which seeds its threads differently from the CPU path
#    NSTAT:   photon number of the statistical comparison (default 4000)
#    SIGMA:   allowed deviation in standard errors (default 4)

MMC=${MMC:-$(dirname "$0")/../bin/mmc}
MMCOPT=${MMCOPT:-"-G -1"}
PYTHON=${PYTHON:-python3}
REFTOL=${REFTOL:-1e-3}
MAXFORK=${MAXFORK:-0}
NSTAT=${NSTAT:-4000}
SIGMA=${SIGMA:-4}
GOLDENTAG=${GOLDENTAG:-}
UPDATE=0

if [ "$1" = "-u" ]; then
    UPDATE=1
    shift
fi

PARAM="$@"
TESTDIR=$(cd "$(dirname "$0")" && pwd)

case "$MMC" in
    */*) MMC=$(cd "$(dirname "$MMC")" && pwd)/$(basename "$MMC") ;;
esac
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
export OMP_NUM_THREADS=1

fail=0

for casefile in "$TESTDIR"/traj_cases/*.json; do
    case=$(basename "$casefile" .json)
    echo "test trajectory of case '$case' ... "

    "$PYTHON" "$TESTDIR/mkmeshcase.py" "$casefile" "$WORKDIR/$case.json"
    (cd "$WORKDIR" && "$MMC" -f "$case.json" $MMCOPT -D M -F jnii -s "$case" $PARAM > "$case.log" 2>&1)

    if [ ! -f "$WORKDIR/${case}_traj.jdt" ]; then
        echo "fail to save trajectory, see log below"
        tail -5 "$WORKDIR/$case.log"
        fail=$((fail + 1))
        continue
    fi

    golden="$TESTDIR/golden/$case$GOLDENTAG.json"

    if [ "$UPDATE" = "1" ] || [ ! -f "$golden" ]; then
        "$PYTHON" "$TESTDIR/trajtest.py" "$WORKDIR/${case}_traj.jdt" --case "$casefile" --wtol "$REFTOL" --save "$golden" --solver mmc || fail=$((fail + 1))
    elif [ "$(uname -s)" != "Linux" ]; then
        # the per-thread RNG seeds come from the C library's srand()/rand(), whose sequence differs
        # between platforms, so the Linux golden paths only apply on Linux; check the invariants here
        "$PYTHON" "$TESTDIR/trajtest.py" "$WORKDIR/${case}_traj.jdt" --case "$casefile" --wtol "$REFTOL" || fail=$((fail + 1))
    else
        "$PYTHON" "$TESTDIR/trajtest.py" "$WORKDIR/${case}_traj.jdt" --case "$casefile" --ref "$golden" --tol "$REFTOL" --wtol "$REFTOL" --max-fork "$MAXFORK" || fail=$((fail + 1))
    fi

    stats="$TESTDIR/golden/${case}_mcxstats.json"

    if [ -f "$stats" ]; then
        echo "compare $NSTAT-photon path statistics with mcx ... "
        "$PYTHON" -c "import json,sys; d=json.load(open(sys.argv[1])); d['Session']['Photons']=int(sys.argv[3]); d['Session']['ID']+='_stat'; json.dump(d,open(sys.argv[2],'w'))" "$WORKDIR/$case.json" "$WORKDIR/${case}_stat.json" "$NSTAT"
        (cd "$WORKDIR" && "$MMC" -f "${case}_stat.json" $MMCOPT -D M -F jnii -s "${case}_stat" $PARAM > "${case}_stat.log" 2>&1)
        "$PYTHON" "$TESTDIR/trajstats.py" "$stats" "$WORKDIR/${case}_stat_traj.jdt" --sigma "$SIGMA" || fail=$((fail + 1))
    fi
done

if [ "$fail" -gt 0 ]; then
    echo "failed $fail trajectory tests"
    exit $fail
fi

echo "passed all trajectory tests!"
