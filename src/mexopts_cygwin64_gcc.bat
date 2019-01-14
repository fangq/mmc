@echo off

set MATLAB=%MATLAB%
set MW_TARGET_ARCH=win64
set PATH=C:\cygwin64\bin;%PATH%

set COMPILER=x86_64-w64-mingw32-g++
set COMPFLAGS=-c -m64 -mwin32 -mdll -Wall -std=c++11 -DMATLAB_MEX_FILE
set OPTIMFLAGS=-DNDEBUG -O2
set DEBUGFLAGS=-g
set NAME_OBJECT=-o

set LINKER=x86_64-w64-mingw32-g++
set LINKFLAGS=-shared -L"%MATLAB%\extern\lib\win64\microsoft" -L"%MATLAB%\bin\win64"
set LINKFLAGSPOST=-lmx -lmex -lmat
set LINKOPTIMFLAGS=-O2
set LINKDEBUGFLAGS=-g
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=-o "%OUTDIR%%MEX_NAME%%MEX_EXT%"