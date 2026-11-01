"""pack_wheel.py <module folder> <pmmc folder> <output folder>

Packs a pmmc wheel for the running Python from the pmmc Python package (pmmc/pmmc)
and the compiled files in <module folder> (_rtmmc*.pyd and the DLLs it loads),
which go at the top of site-packages, as _pmmc does.  build_mmc_trinity.sh builds
_rtmmc with MSVC, which setup.py's CMake build does not use on Windows.
"""

import base64
import hashlib
import os
import re
import sys
import sysconfig
import zipfile

moddir, pkgdir, outdir = sys.argv[1:4]

with open(os.path.join(pkgdir, "setup.py")) as f:
    version = re.search(r'version="([^"]+)"', f.read()).group(1)

tag = "cp%d%d" % sys.version_info[:2]
plat = sysconfig.get_platform().replace("-", "_").replace(".", "_")
wheel = "pmmc-%s-%s-%s-%s.whl" % (version, tag, tag, plat)
distinfo = "pmmc-%s.dist-info" % version

files = []  # (name in the wheel, bytes)
for name in sorted(os.listdir(os.path.join(pkgdir, "pmmc"))):
    if name.endswith(".py"):
        with open(os.path.join(pkgdir, "pmmc", name), "rb") as f:
            files.append(("pmmc/" + name, f.read()))
for name in sorted(os.listdir(moddir)):
    with open(os.path.join(moddir, name), "rb") as f:
        files.append((name, f.read()))

files.append(
    (
        distinfo + "/METADATA",
        (
            "Metadata-Version: 2.1\nName: pmmc\nVersion: %s\n"
            "Summary: Python bindings for Mesh-based Monte Carlo (MMC), with the CUDA and OptiX backends (mmc-trinity)\n"
            "Home-page: https://github.com/fangq/mmc\nAuthor: Qianqian Fang\nLicense: GPLv3+\n"
            "Requires-Python: >=3.6\nRequires-Dist: numpy\n" % version
        ).encode(),
    )
)
files.append(
    (
        distinfo + "/WHEEL",
        (
            "Wheel-Version: 1.0\nGenerator: pmmc pack_wheel.py\nRoot-Is-Purelib: false\nTag: %s-%s-%s\n"
            % (tag, tag, plat)
        ).encode(),
    )
)
files.append((distinfo + "/top_level.txt", b"_rtmmc\npmmc\n"))


def digest(data):
    h = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=")
    return "sha256=" + h.decode()


record = "".join("%s,%s,%d\n" % (n, digest(d), len(d)) for n, d in files)
record += distinfo + "/RECORD,,\n"
files.append((distinfo + "/RECORD", record.encode()))

os.makedirs(outdir, exist_ok=True)
with zipfile.ZipFile(os.path.join(outdir, wheel), "w", zipfile.ZIP_DEFLATED) as z:
    for name, data in files:
        z.writestr(name, data)
print("packed", wheel)
