function buildmmc(varargin)
%
% Format:
%    buildmmc or buildmmc('option1',value1,'option2',value2,...)
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
%      'filelist': a user-defined list of source file names
%
% Dependency (Windows only):
%  1.If you have MATLAB R2017b or later, you may skip this step.
%    To compile mmclab in MATLAB R2017a or earlier on Windows, you must 
%    pre-install the MATLAB support for MinGW-w64 compiler 
%    https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler
%
%    Note: it appears that installing the above Add On is no longer working
%    and may give an error at the download stage. In this case, you should
%    install MSYS2 from https://www.msys2.org/. Once you install MSYS2,
%    run MSYS2.0 MinGW 64bit from Start menu, in the popup terminal window,
%    type
%
%       pacman -Syu
%       pacman -S base-devel gcc git mingw-w64-x86_64-opencl-headers
%
%    Then, start MATLAB, and in the command window, run
%
%       setenv('MW_MINGW64_LOC','C:\msys64\usr');
%
%  2.After installation of MATLAB MinGW support, you must type 
%    "mex -setup C" in MATLAB and select "MinGW64 Compiler (C)". 
%  3.Once you select the MingW C compiler, you should run "mex -setup C++"
%    again in MATLAB and select "MinGW64 Compiler (C++)" to compile C++.
%  4.File C:\Windows\System32\OpenCL.dll must exist. You can obtain this
%    file by installing your graphics driver or install CUDA/AMD GPU SDK
%    and copy OpenCL.dll to the C:\Windows\System32 folder.
%
% This function is part of Mesh-based Monte Carlo (MMC) URL: http://mcx.space
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

cd(fileparts(which(mfilename)));
if(nargin==1 && strcmp(varargin{1},'clean'))
    if(~isempty(dir('*.o')))
        delete('*.o'); 
    end
    return;
end
opt=struct(varargin{:});
pname='mmc';

cflags='-c -Wall -g -DMCX_EMBED_CL -fno-strict-aliasing -m64 -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1 -O3 -fopenmp  -DUSE_OS_TIMER -DUSE_OPENCL';

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
    linkflags='$LINKLIBS -fopenmp -lstdc++ -static';
    cflags=[cflags ' -I./mingw64/include -I"$MINGWROOT/opt/include"'];
    linkflags=[linkflags ' ''C:\Windows\System32\OpenCL.dll'' '];
    linkvar='LINKLIBS';
else
    linkflags='\$CLIBS -fopenmp -static-libgcc -static-libstdc++';
    cflags=[cflags ' -fPIC '];
    linkflags=[linkflags ' -lOpenCL '];
    linkvar='CLIBS';
end
if(~exist('OCTAVE_VERSION','builtin'))
    for i=1:length(filelist)
        flag='CFLAGS';
        if(regexp(filelist{i},'\.[Cc][Pp][Pp]$'))
            flag='CXXFLAGS';
        end
        disp(sprintf('mex OBJEXT=.o %s=''%s'' -c ''%s'' ',flag,cflags,filelist{i}));
        eval(sprintf('mex OBJEXT=.o %s=''%s'' -c ''%s'' ',flag,cflags,filelist{i}));
    end
    if(isfield(opt,'lib'))
        linkflags=[linkflags ' ' opt.lib];
    end
    fn=dir('*.o');
    fprintf(1,  'mex %s -output %s -outdir ../%slab %s=''%s''\n',strjoin({fn.name}),pname,pname,linkvar,linkflags);
    eval(sprintf('mex %s -output %s -outdir ../%slab %s=''%s'' ',strjoin({fn.name}),pname,pname,linkvar,linkflags));
else
    linkflags=regexprep(linkflags,['[\\]*\$' linkvar],'');
    for i=1:length(filelist)
        cmd=sprintf('mex %s -c ''%s'' ',cflags,filelist{i});
        fprintf(stdout,'%s\n',cmd);
	fflush(stdout);
        eval(cmd);
    end
    fn=dir('*.o');
    cmd=sprintf('mex %s -o ../%slab/%s %s ',strjoin({fn.name}),pname,pname,linkflags);
    fprintf(1,'%s\n',cmd);
    fflush(1);
    eval(cmd);
end
