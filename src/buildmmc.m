function buildmmc(varargin)
%
% Format:
%    buildmex or buildmex('option1',value1,'option2',value2,...)
%
% Compiling script for mmclab mex file in MATLAB and GNU Octave. 
% If compiled successfully, the output mex file can be found in the 
% mmc/mmclab folder (or ../mmclab using relative path from mmc/src)
%
% Author: Qianqian Fang <q.fang at neu.edu>
%
% Input:
%    options: without any option, this script compiles mmc.mex* using
%    default settings. Supported options include
%      'include': a string made of sequences of ' -I"/folder/path" ' that 
%            can be included for compilation (format similar to the -I
%            option for gcc)
%      'lib': a string made of sequences of ' -L"/folder/path" ' and '
%           -llibrary' that can be added for linking (format similar to -L
%           and -l flags for gcc)
%
% Dependency (Windows only):
%  1.To compile mmclab in MATLAB R2017a or earlier on Windows, you must 
%    pre-install the MATLAB support for MinGW-w64 compiler 
%    https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler
%  2.After installation of MATLAB MinGW support, you must type "mex -setup"
%    in MATLAB and select "MinGW64 Compiler (C)". 
%  3.File C:\Windows\System32\OpenCL.dll must exist. You can obtain this
%    file by installing your graphics driver or install CUDA/AMD GPU SDK
%    and copy OpenCL.dll to the C:\Windows\System32 folder.
%
% This function is part of Mesh-based Monte Carlo (MMC) URL: http://mcx.space
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

cd(fileparts(which(mfilename)));
if(~isempty(dir('*.o')))
    delete('*.o'); 
end
opt=struct(varargin{:});
pname='mmc';

cflags='-c -Wall -g -DMCX_EMBED_CL -fno-strict-aliasing -m64 -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1 -O3 -fopenmp  -DUSE_OS_TIMER -DUSE_OPENCL';
linkflags='$LINKLIBS -fopenmp -lstdc++ -static';

filelist={'xorshift128p_rand.c','simpmesh.c','tettracing.c',...
    'mcx_utils.c','tictoc.c','cjson/cJSON.c','mmc_host.c',...
    'highordermesh.cpp','mmc_cl_utils.c','mmc_cl_host.c','mmclab.cpp'};
if(isfield(opt,'filelist'))
    filelist=opt.filelist;
end
if(isfield(opt,'include'))
    cflags=[cflags ' ' opt.include];
end
if(ispc)
    cflags=[cflags ' -I./mingw64/include -I"$MINGWROOT/opt/include"'];
    linkflags=[linkflags ' "C:\Windows\System32\OpenCL.dll"'];
end
if(~exist('OCTAVE_VERSION','builtin'))
    for i=1:length(filelist)
        flag='CFLAGS';
        if(regexp(filelist{i},'\.[Cc][Pp][Pp]$'))
            flag='CXXFLAGS';
        end
        eval(sprintf('mex OBJEXT=.o %s=''%s'' -c ''%s'' ',flag,cflags,filelist{i}));
    end
    if(isfield(opt,'lib'))
        linkflags=[linkflags ' ' opt.lib];
    end
    eval(sprintf('mex *.o -output %s -outdir ../%slab LINKLIBS=''%s'' ',pname,pname,linkflags));
else
    for i=1:length(filelist)
        cmd=sprintf('mex %s -c ''%s'' ',cflags,filelist{i});
        disp(cmd);
        eval(cmd);
    end
    eval(sprintf('mex *.o -o ../%slab/%s LINKLIBS=''%s'' ',pname,pname,linkflags));
end
