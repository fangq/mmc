%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Validate DMMC in a heterogeneous domain: sphere in a cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation

% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')

% create a surface mesh for a 10mm radius sphere
[nsph,fsph]=meshasphere([30 30 30],10,1.0);
[nsph,fsph]=removeisolatednode(nsph,fsph);

%% create DMMC and MMC simulation parameters
clear cfg

cfg.nphoton=1e7;
cfg.seed=1648335518;
cfg.srcpos=[30,30.001,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;
cfg.prop=[0 0 1 1;0.002 1.0 0.01 1.37;0.050 0 1 1.37];
cfg.debuglevel='TP';
cfg.isreflect=1;

cfgs=[cfg cfg];

%% tetrahedral mesh generation
% generate a coarse volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 20

ISO2MESH_SESSION='dmmc_';

[nbox,ebox]=meshgrid6(0:60:60,0:60:60,0:60:60);
fbox=volface(ebox);
[no,fc]=mergemesh(nbox,fbox,nsph,fsph);

ISO2MESH_TETGENOPT='-Y -A'
[cfgs(1).node,cfgs(1).elem]=surf2mesh(no,fc,[0 0 0],[60.1 60.1 60.1],1,100,[1,1,1;30 30 30]);
cfgs(1).method='grid';

clear ISO2MESH_TETGENOPT
ISO2MESH_SESSION='mmc_';

[no,el]=meshasphere([30 30 30],10,1.0);
nodesize=[1*ones(size(no,1),1) ; 1; 1];
srcpos=[30. 30.001 0.];
fixednodes=[30.,30.,0.05; 30 30 30];
nfull=[no;fixednodes];
[cfgs(2).node,cfgs(2).elem,face2]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,2,[30 30 30],[],[1 1 1 1 1 1 1 1]);
[cfgs(2).node,cfgs(2).elem]=sortmesh(srcpos,cfgs(2).node,cfgs(2).elem,1:4);
cfgs(2).elem(:,5)=cfgs(2).elem(:,5)+1;
cfgs(2).method='elem';

%% running simulation with 1x1x1 mm input volume (60x60x60 grid)
phimmc=mmclab(cfgs);

%% create MCX simulation parameters
clear cfg;

% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E');

cfg.nphoton=3e7;

% define a 1cm radius sphere within a 6x6x6 cm box
dim=60;
[xi,yi,zi]=meshgrid(0.5:dim-0.5,0.5:dim-0.5,0.5:dim-0.5);
dist=(xi-30).^2+(yi-30).^2+(zi-30).^2;
cfg.vol=ones(size(xi));
cfg.vol(dist<100)=2;
cfg.vol=uint8(cfg.vol);

% define the source position, same as DMMC/MMC
cfg.srcpos=[30,30.001,0];
cfg.srcdir=[0 0 1];
cfg.issrcfrom0=1;

% format: [mua(1/mm) mus(1/mm) g n]
cfg.prop=[0 0 1 1;0.002 1.0 0.01 1.37;0.050 0 1 1.37];

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;
  
% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

% enable reflection at exterior boundaries, same as DMMC/MMC
cfg.isreflect=1; 

%% running simulation with 1x1x1 mm input volume (60x60x60 grid)
fprintf('running MCX simulation ... this takes about 6 seconds on a GTX 470\n');
tic;
phimcx=mcxlab(cfg);
toc;

%% visualization: compare CW fluence maps produced by DMMC/MMC against that by MCX using contour lines
clines=-7:0.25:-1;
[xx,zz]=meshgrid(0.5:59.5,0.5:59.5);
node2=cfgs(2).node;
elem2=cfgs(2).elem;
flux_mmc=phimmc(2).data;
CWfluence_mmc=sum(flux_mmc,2)*cfg.tstep;
[cutpos,cutvalue,facedata]=qmeshcut(elem2(:,1:4),node2,CWfluence_mmc,[0 30.5 0; 0 30.5 1; 1 30.5 0]);

CWfluence_mmc_grid=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xx,zz);
CWfluence_dmmc=sum(phimmc(1).data,4)*cfg.tstep;
CWfluence_mcx=sum(phimcx.data,4)*cfg.tstep;

figure;
contourf(xx,zz,log10(abs(squeeze(CWfluence_mmc_grid))),clines);
hold on
contour(xx,zz,log10(abs(squeeze(CWfluence_dmmc(1:end-1,31,1:end-1))')),clines,'linestyle','--','color','w','linewidth',1.5);
contour(xx,zz,log10(abs(squeeze(CWfluence_mcx(:,31,:))')),clines,'linestyle','--','color',[0.9100    0.4100    0.1700],'linewidth',1.5);

[xcirc,ycirc] = cylinder(10,200);
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal;
colormap;
set(gca,'fontsize',20);
lg=legend('MMC','DMMC','MCX');
set(lg,'color','none');
set(lg,'box','off');
