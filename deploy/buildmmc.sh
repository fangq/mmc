#!/bin/bash

###############################################################################
#
#  MMC Nightly-Build Script
#
#  by Qianqian Fang <q.fang at neu.edu>
#
#  Format:
#     ./buildmmc.sh <releasetag> <branchname>
#                   releasetag defaults to "nightly" if not given
#                   branchname defaults to the default main branch
#
#  Dependency:
#   - To compile mmc binary, mmclab for octave
#
#     for Linux:
#         sudo apt-get install gcc liboctave-dev upx-ucl vim-common ocl-icd-opencl-dev
#
#   - To compile mmclab for MATLAB, one must install MATLAB first, also search
#     and replace R20xx in this script to match your system's MATLAB version
#   - For Windows, first install Cygwin64, and install x86_64-w64-mingw32-gcc/g++
#     or install MSYS2 with mingw64 based gcc compilers
#
###############################################################################


## setting up environment

BUILD='nightly'
if [ ! -z "$1" ]; then
	BUILD=$1
fi
DATE=$(date +'%Y%m%d')
BUILDROOT=~/space/autobuild/$BUILD/mmc
OSID=$(uname -s)
MACHINE=$(uname -m)

if [ "$OSID" == "Linux" ]; then
	OS=linux
	source ~/.bashrc
elif [ "$OSID" == "Darwin" ]; then
	OS=macos
	source ~/.bash_profile
elif [[ $OSID == CYGWIN* ]] || [[ $OSID == MINGW* ]] || [[ $OSID == MSYS* ]]; then
	OS=win
fi

TAG=${OS}-${MACHINE}-${BUILD}

if [ "$BUILD" == "nightly" ]; then
	TAG=${OS}-${MACHINE}-${BUILD}build
fi

## setting up upload server (blank if no need to upload)

SERVER=
REMOTEPATH=

## checking out latest github code

mkdir -p $BUILDROOT
cd $BUILDROOT

rm -rf mmc
git clone https://github.com/fangq/mmc.git

## automatically update revision/version number

cat <<EOF >>mmc/.git/config
[filter "rcs-keywords"]
        clean  = .git_filters/rcs-keywords.clean
        smudge = .git_filters/rcs-keywords.smudge %f
EOF

cd mmc
if [ ! -z "$2" ]; then
	git checkout $2
fi
rm -rf *
git checkout *
rm -rf .git

## zip and upload source code package

cd ..
zip -FSr $BUILDROOT/mmc-src-${BUILD}.zip mmc
if [ "$OS" == "linux" ] && [ ! -z "$SERVER" ]; then
	scp $BUILDROOT/mmc-src-${BUILD}.zip $SERVER:$REMOTEPATH/src/
fi

## build matlab mex file

cd mmc/src

rm -rf ../mmclab/AUTO_BUILD_*
make clean
if [ "$OS" == "win" ]; then
	make mex &> ../mmclab/AUTO_BUILD_${DATE}.log
elif [ "$OS" == "osx" ]; then
	make mex &> ../mmclab/AUTO_BUILD_${DATE}.log
else
	make mex MEXLINKOPT="-static-libstdc++ -static-libgcc -fopenmp" &>../mmclab/AUTO_BUILD_${DATE}.log
fi

## build octave mex file

make clean
if [ "$OS" == "macos" ]; then
	make oct USEROCTOPT="CXXFLAGS='-pipe -Os -arch x86_64' DL_LD=g++ DL_LDFLAGS='-fopenmp -static-libgcc -static-libstdc++'" >>../mmclab/AUTO_BUILD_${DATE}.log 2>&1
elif [ "$OS" == "win" ]; then
	OLDPATH="$PATH"
	export PATH="C:\Octave\Octave-8.2.1\mingw64\bin":$PATH
	make oct CC=gcc MEXLINKOPT='"C:\tools\msys64\mingw64\lib\gcc\x86_64-w64-mingw32\10.3.0\libgomp.a" "C:\cygwin64\usr\x86_64-w64-mingw32\sys-root\mingw\lib\libwinpthread.a" -static-libgcc -static-libstdc++' >>../mmclab/AUTO_BUILD_${DATE}.log 2>&1
	export PATH="$OLDPATH"
else
	make oct MEXLINKOPT="-static-libgcc -static-libstdc++ -Wl,-Bstatic -lm -lpthread -Wl,-Bdynamic" >>../mmclab/AUTO_BUILD_${DATE}.log 2>&1
fi

## test mex file dependencies

mexfile=(../mmclab/mmc.mex*)

if [ -f "$mexfile" ]; then
	if [ "$OS" == "macos" ]; then
		otool -L ../mmclab/mmc.mex* >>../mmclab/AUTO_BUILD_${DATE}.log
	elif [ "$OS" == "win" ]; then
		objdump -p ../mmclab/mmc.mex* | grep -H "DLL Name:" >>../mmclab/AUTO_BUILD_${DATE}.log
	else
		ldd ../mmclab/mmc.mex* >>../mmclab/AUTO_BUILD_${DATE}.log
	fi
	echo "Build Octave MMCLAB Successfully" >>../mmclab/AUTO_BUILD_${DATE}.log
else
	echo "Build Octave MMCLAB Failed" >>../mmclab/AUTO_BUILD_${DATE}.log
fi

## compress mex files with upx

upx -9 ../mmclab/mmc.mex* || true

## zip and upload mex package

if [ "$BUILD" != "nightly" ]; then
	rm -rf ../mmclab/AUTO_BUILD_${DATE}.log
fi

cp $BUILDROOT/dlls/*.dll ../mmclab
cd ..
zip -FSr $BUILDROOT/mmclab-${TAG}.zip mmclab
cd src
[ ! -z "$SERVER" ] && scp $BUILDROOT/mmclab-${TAG}.zip $SERVER:$REMOTEPATH/${OS}64/

## compile standalone binary/executable

make clean

if [ "$OS" == "macos" ]; then
	make &>$BUILDROOT/mmc_buildlog_${DATE}.log
elif [ "$OS" == "win" ]; then
	make &>$BUILDROOT/mmc_buildlog_${DATE}.log
else
	make AR=c++ EXTRALIB="-static-libstdc++ -static-libgcc -lOpenCL -lm" &>$BUILDROOT/mmc_buildlog_${DATE}.log
fi

## test binary dependencies

if [ -f "../bin/mmc" ]; then
	if [ "$OS" == "macos" ]; then
		otool -L ../bin/mmc >>$BUILDROOT/mmc_buildlog_${DATE}.log
	elif [ "$OS" == "win" ]; then
		objdump -p ../bin/mmc.exe | grep "DLL Name:" >>$BUILDROOT/mmc_buildlog_${DATE}.log
	else
		ldd ../bin/mmc >>$BUILDROOT/mmc_buildlog_${DATE}.log
	fi
	echo "Build MMC Binary Successfully" >>$BUILDROOT/mmc_buildlog_${DATE}.log
else
	echo "Build MMC Binary Failed" >>$BUILDROOT/mmc_buildlog_${DATE}.log
	exit 1
fi

## compress binary with upx

upx -9 ../bin/mmc* || true

## zip and upload binary package

cd ../
rm -rf .git* .travis* mmclab webmmc commons win32 pmcxcl deploy
#cp $BUILDROOT/dlls/*.dll bin
rm -rf src .git_filters .gitattributes
cd ../

mv $BUILDROOT/mmc_buildlog_${DATE}.log mmc/AUTO_BUILD_${DATE}.log

if [ "$BUILD" != "nightly" ]; then
	rm -rf mmc/AUTO_BUILD_${DATE}.log
fi

if [ "$OS" == "win" ]; then
	zip -FSr mmc-${TAG}.zip mmc
else
	zip -FSry mmc-${TAG}.zip mmc
fi

#mv mmc-${TAG}.zip $BUILDROOT

cd $BUILDROOT

[ ! -z "$SERVER" ] && scp mmc-${TAG}.zip $SERVER:$REMOTEPATH/${OS}64/
