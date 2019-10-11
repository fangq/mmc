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
cfg.srcpos=[30.5,30.5,0];
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
cfgs(1).method='grid';

% generate a refined volumetric mesh from triangular surfaces of 
% three spheres with an additional bounding box for MMC
clear ISO2MESH_TETGENOPT
ISO2MESH_SESSION='dmmc2_';

[no,fc]=mergemesh(nsph_10,fsph_10,nsph_23,fsph_23,nsph_25,fsph_25);
[no,fc]=removeisolatednode(no,fc);
srcpos=[30.5 30.5 0.];
[cfgs(2).node,cfgs(2).elem,face2]=surf2mesh(no,fc,[0 0 0],[60 60 60],1,0.25,[1 1 1 0.1;30 30 6 0.1;30 30 17 0.1;30 30 30 0.1],[],[1 1 1 1 1 1 1 1]);
[cfgs(2).node,cfgs(2).elem]=sortmesh(srcpos,cfgs(2).node,cfgs(2).elem,1:4);
cfgs(2).method='elem';

% % 3D view of the cross-section of the mesh: 'y>30'
% figure;plotmesh(cfgs(1).node,cfgs(1).elem,'y>30');view(3);
% figure;plotmesh(cfgs(2).node,cfgs(2).elem,'y>30');view(3);

phimmc=mmclab(cfgs);

%% mcx
clear cfg;

% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E');

cfg.nphoton=1e8;

% define three spheres with radius=10 mm, 23 mm and 25 mm within a 60x60x60 mm box
dim=60;
[xi,yi,zi]=meshgrid(0.5:(dim-0.5),0.5:(dim-0.5),0.5:(dim-0.5));
dist=(xi-30).^2+(yi-30).^2+(zi-30).^2;
cfg.vol=ones(size(xi));
cfg.vol(dist<625)=2;
cfg.vol(dist<529)=3;
cfg.vol(dist<100)=4;
cfg.vol=uint8(cfg.vol);

% define the source position
cfg.srcpos=[30.5,30.5,0];
cfg.srcdir=[0 0 1];
cfg.issrcfrom0=1;

% format: [mua(1/mm) mus(1/mm) g n]
cfg.prop=[0,0,1,1;0.02 7.0 0.89 1.37;0.004 0.009 0.89 1.37;0.02 9.0 0.89 1.37;0.05 0.0 1.0 1.37];

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;
  
% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

cfg.isreflect=1; % disable reflection at exterior boundaries
cfg.unitinmm=1; % resolution

%% running simulation with 1x1x1 mm input volume (60x60x60 grid)
fprintf('running simulation ... this takes about 6 seconds on a GTX 470\n');
tic;
phimcx=mcxlab(cfg);
toc;

%% data visualization
%mcx
phi_mcx=sum(phimcx.data,4);
%mmc
[xx,yy]=meshgrid(0.5:59.5,0.5:59.5);
node1=cfgs(2).node;
elem1=cfgs(2).elem;
mesh0=phimmc(2).data;
s1=sum(mesh0,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem1(:,1:4),node1,s1,[0 30.5 0; 0 30.5 1; 1 30.5 0]);
phi_mmc=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xx,yy);
%dmmc
phi_dmmc=sum(phimmc(1).data,4);
phi_dmmc=phi_dmmc(1:60,1:60,1:60); 

figure;
clines=[-10:0.5:10];
contourf(log10(abs(squeeze(phi_mmc))),clines,'linewidth',1.5);
hold on
contour(log10(abs(squeeze(phi_dmmc(:,31,:))')),clines,'linestyle','--','color','w','linewidth',1.5);
contour(log10(abs(squeeze(phi_mcx(:,31,:))')),clines,'linestyle','--','color',[0.9100    0.4100    0.1700],'linewidth',1.5);
colorbar('EastOutside');

%plot a dashed circle with radius of 10
[xcirc,ycirc] = cylinder(10,200);
xcirc=xcirc(1,:)+30.5;
ycirc=ycirc(1,:)+30.5;
plot(xcirc,ycirc,'k--','linewidth',2);

%plot a dashed circle with radius of 23
[xcirc,ycirc] = cylinder(23,200);
xcirc=xcirc(1,:)+30.5;
ycirc=ycirc(1,:)+30.5;
plot(xcirc,ycirc,'k--','linewidth',2);

%plot a dashed circle with radius of 25
[xcirc,ycirc] = cylinder(25,200);
xcirc=xcirc(1,:)+30.5;
ycirc=ycirc(1,:)+30.5;
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal;
colormap;
lg=legend('MMC','DMMC','MCX','Location','northoutside','orientation','horizontal');
set(lg,'Color',[0.5 0.5 0.5]);
set(lg,'box','on');
set(gca,'fontsize',18);