#!/bin/sh
#
# testmmc.sh - functional tests for MMC
#
# usage: ./testmmc.sh [extra mmc parameters]
#
# runs the built-in benchmarks on the SSE/OpenMP CPU path (-G -1) and checks
# the absorbed fraction against the expected values, then runs the
# step-by-step trajectory tests in testtraj.sh. The GPU backends compiled
# into the binary (opencl; cuda with "make cuda"; optix with "make optix")
# and usable on device -G 1 are then tested as well: absorbed fraction of the
# dual-grid benchmarks for opencl/cuda, step-by-step trajectories of cuda
# against the opencl golden data, and for optix (which does not report the
# absorbed fraction or save trajectories) a smoke run plus a comparison of the
# planar-source fluence with opencl.
#
# environment variables:
#    MMC:  path to the mmc binary (default ../bin/mmc)
#    BACKENDS: GPU backends to test (default: probe opencl, cuda and optix)
#    PYTHON: python with numpy, used by the optix planar test (default python3)

MMC=${MMC:-$(dirname "$0")/../bin/mmc}
PYTHON=${PYTHON:-python3}
PARAM=$@
TESTDIR=$(cd "$(dirname "$0")" && pwd)

case "$MMC" in
    */*) MMC=$(cd "$(dirname "$MMC")" && pwd)/$(basename "$MMC") ;;
esac

fail=0

echo "test binary file ... "
if [ ! -f "$MMC" ]; then echo "mmc binary $MMC does not exist"; exit 1; else echo "ok"; fi

echo "test version number ... "
temp=$("$MMC" --version 2>&1 | grep -E 'v20[0-9]+')
if [ -z "$temp" ]; then echo "fail to print version"; fail=$((fail + 1)); else echo "ok"; fi

echo "test help info ... "
temp=$("$MMC" 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -c -E '^ +-')
if [ "$temp" -lt 40 ]; then echo "fail to print help info"; fail=$((fail + 1)); else echo "ok"; fi

echo "test listing devices ... "
temp=$("$MMC" -L 2>&1 | grep -i -E 'platform|device|SSE')
if [ -z "$temp" ]; then echo "fail to list devices"; fail=$((fail + 1)); else echo "ok"; fi

echo "test exporting json input ... "
temp=$("$MMC" --bench dmmc-cube60 --dumpjson 1 2>&1 | grep -E '"MeshNode"')
if [ -z "$temp" ]; then echo "fail to export json input"; fail=$((fail + 1)); else echo "ok"; fi

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# benchmark name and the expected integer part(s) of the absorbed fraction
for bench in "dmmc-cube60:1[78]" "dmmc-cube60b:2[67]" "edgeimmc:30" "nodeimmc:1[45]" "faceimmc:94"; do
    name=${bench%%:*}
    expect=${bench#*:}
    echo "test absorbed fraction of benchmark $name ... "
    temp=$(cd "$WORKDIR" && "$MMC" --bench "$name" -G -1 -n 1e5 -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E "absorbed: *$expect\.[0-9]+%")
    if [ -z "$temp" ]; then echo "fail to match the absorbed fraction of $name"; fail=$((fail + 1)); else echo "ok"; fi
done

# the ROI medium (mua=23.05) of edgeimmc shares g=0.9 with the background;
# setting only its g to 0 must change the result of a fixed-seed run
echo "test implicit ROI scattering uses the anisotropy of the ROI medium ... "
a=""
for roig in keep 0; do
    if [ "$roig" = keep ]; then filter=cat; else filter='sed -E /"mua":[[:space:]]*23\.05/{n;n;s/"g":[[:space:]]*[0-9.]+/"g":0/;}'; fi
    a="$a $(cd "$WORKDIR" && "$MMC" --bench edgeimmc -G -1 --dumpjson - 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | $filter |
            OMP_NUM_THREADS=1 "$MMC" -f - -G -1 -n 1e4 -E 12345 -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E 'absorbed: *[0-9.]+%' | tr -d ' ')"
done
set -- $a
if [ $# -ne 2 ] || [ "$1" = "$2" ]; then echo "fail: ROI g has no effect ($a)"; fail=$((fail + 1)); else echo "ok"; fi

echo "test photon replay ... "
temp=$(cd "$WORKDIR" && "$MMC" --bench dmmc-cube60b -G -1 -n 1e5 -S 0 -q 1 -s replaytest -D T $PARAM > /dev/null 2>&1 &&
       "$MMC" --bench dmmc-cube60b -G -1 -E replaytest_detp.jdt -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E 'simulated.*photons|absorbed: *[0-9.]+%')
if [ -z "$temp" ]; then echo "fail to replay detected photons"; fail=$((fail + 1)); else echo "ok"; fi

echo "test step-by-step photon trajectories ... "
if ! MMC="$MMC" sh "$TESTDIR/testtraj.sh" > "$WORKDIR/traj.log" 2>&1; then
    cat "$WORKDIR/traj.log"
    fail=$((fail + 1))
else
    echo "ok"
fi

if [ -z "${BACKENDS+x}" ]; then
    BACKENDS=""
    # a backend missing from the binary falls back to another one, so look for its own log line
    for b in opencl cuda optix; do
        case $b in
            opencl) mark='OpenCL' ;;
            cuda)   mark='compiled by nvcc' ;;
            optix)  mark='optix init complete' ;;
        esac
        log=$(cd "$WORKDIR" && "$MMC" --bench dmmc-cube60 -c $b -G 1 -n 1000 -S 0 -D T 2>&1)
        if echo "$log" | grep -q 'simulated 1000 photons' && echo "$log" | grep -q "$mark"; then
            BACKENDS="$BACKENDS $b"
        fi
    done
fi

for b in $BACKENDS; do
    if [ "$b" = "optix" ]; then
        echo "test running benchmark dmmc-cube60 with optix ... "
        temp=$(cd "$WORKDIR" && "$MMC" --bench dmmc-cube60 -c optix -G 1 -n 1e5 -S 0 -D T $PARAM 2>&1 | grep -o -E 'simulated 100000 photons')
        if [ -z "$temp" ]; then echo "fail to run optix"; fail=$((fail + 1)); else echo "ok"; fi

        # a planar source launched in the air above the mesh must not deposit
        # fluence in the gap before it enters the tissue
        case " $BACKENDS " in
            *" opencl "*)
                echo "test planar source fluence of optix against opencl ... "
                cp "$TESTDIR"/../examples/planar/planar.inp "$TESTDIR"/../examples/planar/*_digimouse.dat "$WORKDIR"/
                for c in optix opencl; do
                    (cd "$WORKDIR" && "$MMC" -f planar.inp -b 1 -n 1e5 -c $c -G 1 -s planar_$c -F bin $PARAM > planar_$c.log 2>&1)
                done
                if ! (cd "$WORKDIR" && "$PYTHON" -c "import numpy as np; a=np.fromfile('planar_optix.bin').sum(); b=np.fromfile('planar_opencl.bin').sum(); print('optix/opencl total fluence ratio %.4f' % (a/b)); raise SystemExit(int(not abs(a/b-1)<0.05))"); then
                    fail=$((fail + 1))
                else
                    echo "ok"
                fi
                ;;
        esac
        continue
    fi

    for bench in "dmmc-cube60:1[78]" "dmmc-cube60b:2[67]"; do
        name=${bench%%:*}
        expect=${bench#*:}
        echo "test absorbed fraction of benchmark $name with $b ... "
        temp=$(cd "$WORKDIR" && "$MMC" --bench "$name" -c $b -G 1 -n 1e5 -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E "absorbed: *$expect\.[0-9]+%")
        if [ -z "$temp" ]; then echo "fail to match the absorbed fraction of $name with $b"; fail=$((fail + 1)); else echo "ok"; fi
    done

    if [ "$b" = "cuda" ]; then
        echo "test step-by-step photon trajectories with cuda against the opencl golden data ... "
        if ! MMC="$MMC" MMCOPT="-G 1 -c cuda" GOLDENTAG=_opencl REFTOL=5e-2 MAXFORK=0.1 sh "$TESTDIR/testtraj.sh" -t 4000 -T 1 > "$WORKDIR/trajcuda.log" 2>&1; then
            cat "$WORKDIR/trajcuda.log"
            fail=$((fail + 1))
        else
            echo "ok"
        fi
    fi
done

if [ "$fail" -gt 0 ]; then
    echo "failed $fail tests"
    exit $fail
fi

echo "passed all tests!"
