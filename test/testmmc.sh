#!/bin/sh
#
# testmmc.sh - functional tests for MMC
#
# usage: ./testmmc.sh [extra mmc parameters]
#
# runs the built-in benchmarks on the SSE/OpenMP CPU path (-G -1) and checks
# the absorbed fraction against the expected values, then runs the
# step-by-step trajectory tests in testtraj.sh
#
# environment variables:
#    MMC:  path to the mmc binary (default ../bin/mmc)

MMC=${MMC:-$(dirname "$0")/../bin/mmc}
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
for bench in "dmmc-cube60:1[78]" "dmmc-cube60b:2[67]" "edgeimmc:30" "nodeimmc:14" "faceimmc:94"; do
    name=${bench%%:*}
    expect=${bench#*:}
    echo "test absorbed fraction of benchmark $name ... "
    temp=$(cd "$WORKDIR" && "$MMC" --bench "$name" -G -1 -n 1e5 -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E "absorbed: *$expect\.[0-9]+%")
    if [ -z "$temp" ]; then echo "fail to match the absorbed fraction of $name"; fail=$((fail + 1)); else echo "ok"; fi
done

echo "test photon replay ... "
temp=$(cd "$WORKDIR" && "$MMC" --bench dmmc-cube60b -G -1 -n 1e5 -S 0 -q 1 -s replaytest -D T $PARAM > /dev/null 2>&1 &&
       "$MMC" --bench dmmc-cube60b -G -1 -E replaytest_detp.jdt -S 0 -D T $PARAM 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -o -E 'simulated.*photons|absorbed: *[0-9.]+%')
if [ -z "$temp" ]; then echo "fail to replay detected photons"; fail=$((fail + 1)); else echo "ok"; fi

echo "test step-by-step photon trajectories ... "
if ! sh "$TESTDIR/testtraj.sh" > "$WORKDIR/traj.log" 2>&1; then
    cat "$WORKDIR/traj.log"
    fail=$((fail + 1))
else
    echo "ok"
fi

if [ "$fail" -gt 0 ]; then
    echo "failed $fail tests"
    exit $fail
fi

echo "passed all tests!"
