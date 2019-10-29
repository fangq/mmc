%%-----------------------------------------------------------------
%% validate DMMC with MMC in a heterogeneous cubic domain
%  (see DMMC paper Fig. 1, and MMCL paper Fig. 1b)
%%-----------------------------------------------------------------

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

cfg.nphoton=1e8;
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
cfgs(1).method='grid';

% generate a refined volumetric mesh from triangular surfaces of 
% three spheres with an additional bounding box for MMC
clear ISO2MESH_TETGENOPT
ISO2MESH_SESSION='dmmc2_';

[no,fc]=mergemesh(nsph_10,fsph_10,nsph_23,fsph_23,nsph_25,fsph_25);
[no,fc]=removeisolatednode(no,fc);
srcpos=cfgs(1).srcpos;
[cfgs(2).node,cfgs(2).elem,face2]=surf2mesh(no,fc,[0 0 0],[60 60 60],1,0.25,[1 1 1 0.1;30 30 6 0.1;30 30 17 0.1;30 30 30 0.1],[],[1 1 1 1 1 1 1 1]);
[cfgs(2).node,cfgs(2).elem]=sortmesh(srcpos,cfgs(2).node,cfgs(2).elem,1:4);
cfgs(2).method='elem';

%%

cfgs(1).gpuid=1;
b2dgpu=mmclab(cfgs(1));
b2dgpu=sum(b2dgpu.data,4);

%%
cfgs(1).gpuid=-1;
b2dcpu=mmclab(cfgs(1));
b2dcpu=sum(b2dcpu.data,4);

%%
cfgs(2).gpuid=1;
b2gpu=mmclab(cfgs(2));
b2gpu=sum(b2gpu.data,2);

%%
cfgs(2).gpuid=-1;
b2cpu=mmclab(cfgs(2));
b2cpu=sum(b2cpu.data,2);

%% 
[xi,yi]=meshgrid(0.5:59.5,0.5:59.5);
[cutpos1,cutvalue1,facedata1]=qmeshcut(cfgs(2).elem(:,1:4),cfgs(2).node,b2cpu,[0 30.5 0; 0 30.5 1; 1 30.5 0]);
[cutpos2,cutvalue2,facedata2]=qmeshcut(cfgs(2).elem(:,1:4),cfgs(2).node,b2gpu,[0 30.5 0; 0 30.5 1; 1 30.5 0]);
b2cpu_interpolated=griddata(cutpos1(:,1),cutpos1(:,3),cutvalue1,xi,yi);
b2gpu_interpolated=griddata(cutpos2(:,1),cutpos2(:,3),cutvalue2,xi,yi);

%% mcxcl - add voxelated mcxcl results

cfg_mcx.nphoton=cfgs(1).nphoton;

% define three spheres with radius=10 mm, 23 mm and 25 mm within a 60x60x60 mm box
dim=60;
[xi,yi,zi]=meshgrid(0.5:(dim-0.5),0.5:(dim-0.5),0.5:(dim-0.5));
dist=(xi-30).^2+(yi-30).^2+(zi-30).^2;
cfg_mcx.vol=ones(size(xi));
cfg_mcx.vol(dist<625)=2;%radius 25
cfg_mcx.vol(dist<529)=3;%radius 23
cfg_mcx.vol(dist<100)=4;%radius 10
cfg_mcx.vol=uint8(cfg_mcx.vol);

% define the source position
cfg_mcx.srcpos=[30,30.1,0];
cfg_mcx.srcdir=[0 0 1];
cfg_mcx.issrcfrom0=1;

% format: [mua(1/mm) mus(1/mm) g n]
cfg_mcx.prop=[0,0,1,1;0.02 7.0 0.89 1.37;0.004 0.009 0.89 1.37;0.02 9.0 0.89 1.37;0.05 0.0 1.0 1.37];

% time-domain simulation parameters
cfg_mcx.tstart=0;
cfg_mcx.tend=5e-9;
cfg_mcx.tstep=5e-10;
  
% GPU thread configuration
cfg_mcx.autopilot=1;
cfg_mcx.gpuid=1;

cfg_mcx.isreflect=1; % disable reflection at exterior boundaries
cfg_mcx.unitinmm=1; % resolution

phimcx=mcxlabcl(cfg_mcx);
b2mcx=sum(phimcx.data,4);

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.5
%%-----------------------------------------------------------------
figure
clines = 0:-0.5:-8;
[xi,yi]=meshgrid(0.5:59.5,0.5:59.5);
srcs=[30.1,30.2,0];
dets=[xi(:) 30.2*ones(size(xi(:))) yi(:)];

hold on
[c h2]=contourf(xi,yi, log10(max(squeeze(b2cpu_interpolated*cfgs(1).tstep),1e-8)), clines, 'k-','linewidth', 2 );
contour(xi,yi,log10(abs(b2gpu_interpolated*cfgs(1).tstep)),clines,'g--','linewidth',2);
contour(xi,yi,log10(abs(squeeze(b2dcpu(1:end-1,31,1:end-1))'*cfgs(1).tstep)),clines,'b--','linewidth',2)
contour(xi,yi,log10(abs(squeeze(b2dgpu(1:end-1,31,1:end-1))'*cfgs(1).tstep)),clines,'r-','linewidth',1)
contour(xi,yi,log10(abs(squeeze(b2mcx(:,31,:))'*cfg_mcx.tstep)),clines,'m--')

%plot a dashed circle with radius of 10
[xcirc,ycirc] = cylinder(sqrt(10^2-0.5^2),200); %since we are ploting at y=30.5,0.5 from center of the sphere.
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;
plot(xcirc,ycirc,'k--','linewidth',2);

%plot a dashed circle with radius of 23
[xcirc,ycirc] = cylinder(sqrt(23^2-0.5^2),200);
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;
plot(xcirc,ycirc,'k--','linewidth',2);

%plot a dashed circle with radius of 25
[xcirc,ycirc] = cylinder(sqrt(25^2-0.5^2),200);
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal  
set(gca,'xlim',[0.5 59.5])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('MMC','MMCL','D-MMC','D-MMCL','MCX-CL')
legend boxoff;
box on;