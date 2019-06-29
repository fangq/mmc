%%-----------------------------------------------------------------
%% validate DMMC with MMC in a heterogeneous cubic domain
%  (see DMMC paper Fig. 1)
%%-----------------------------------------------------------------
%
% In this example, we validate the DMMC algorithm with MMC using a 
% heterogeneous cubic domain. We also compare the DMMC/MMC solution with
% that of voxel-based MC (MCX). This simulation generates the results
% shown in Fig. 1(c-e) in the paper.
% 
%%-----------------------------------------------------------------

% preparation

% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

% 2. you need to add the path to MMC matlab folder
% addpath('../../matlab');

% create a surface mesh for a 10 mm radius sphere
[nsph_10,fsph_10]=meshasphere([30 30 30],10,1.0);
[nsph_10,fsph_10]=removeisolatednode(nsph_10,fsph_10);

% create a surface mesh for a 23 mm radius sphere
[nsph_23,fsph_23]=meshasphere([30 30 30],23,2.0);
[nsph_23,fsph_23]=removeisolatednode(nsph_23,fsph_23);

% create a surface mesh for a 25 mm radius sphere
[nsph_25,fsph_25]=meshasphere([30 30 30],25,2.0);
[nsph_25,fsph_25]=removeisolatednode(nsph_25,fsph_25);

%%-----------------------------------------------------------------
%% create simulation parameters
%%-----------------------------------------------------------------
clear cfg

cfg.nphoton=3e6;
cfg.seed=1648335518;
cfg.srcpos=[30,30.1,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;
cfg.prop=[0,0,1,1;0.02 7.0 0.89 1.37;0.004 0.009 0.89 1.37;0.02 9.0 0.89 1.37;0.05 0.0 1.0 1.37];
cfg.debuglevel='TP';
cfg.isreflect=1;

cfgs=[cfg cfg];

%%-----------------------------------------------------------------
%% tetrahedral mesh generation
%%-----------------------------------------------------------------

% generate a coarse CDT volumetric mesh from triangular surfaces of 
% three spheres with an additional bounding box for DMMC

ISO2MESH_SESSION='dmmc1_';

[nbox,ebox]=meshgrid6(0:60:60,0:60:60,0:60:60);
fbox=volface(ebox);
[no,fc]=mergemesh(nsph_10,fsph_10,nsph_23,fsph_23,nsph_25,fsph_25,nbox,fbox);
[no,fc]=removeisolatednode(no,fc);

ISO2MESH_TETGENOPT='-A -q -Y'
[cfgs(1).node,cfgs(1).elem]=surf2mesh(no,fc,[0 0 0],[60.1 60.1 60.1],1,100,[1 1 1;30 30 6;30 30 15;30 30 30]);%thin layer
% [cfgs(1).node,cfgs(1).elem]=surf2mesh(no,fc,[0 0 0],[60.1 60.1 60.1],1,100,[1 1 1;30 30 6;30 30 17;30 30 30]);%thick layer
cfgs(1).method='g';

% generate a refined volumetric mesh from triangular surfaces of 
% three spheres with an additional bounding box for MMC
clear ISO2MESH_TETGENOPT
ISO2MESH_SESSION='dmmc2_';

[no,fc]=mergemesh(nsph_10,fsph_10,nsph_23,fsph_23,nsph_25,fsph_25);
[no,fc]=removeisolatednode(no,fc);
srcpos=cfgs(1).srcpos;
[cfgs(2).node,cfgs(2).elem,face2]=surf2mesh(no,fc,[0 0 0],[60 60 60],1,0.25,[1 1 1 0.1;30 30 6 0.1;30 30 17 0.1;30 30 30 0.1],[],[1 1 1 1 1 1 1 1]);
[cfgs(2).node,cfgs(2).elem]=sortmesh(srcpos,cfgs(2).node,cfgs(2).elem,1:4);
cfgs(2).method='s';

mmc2json(cfgs(1),'dmmc_sphshells.json');
mmc2json(cfgs(2),'sphshells.json');
