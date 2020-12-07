%%-----------------------------------------------------------------
%% validate MMC with a homogeneous cubic domain
%  (see MMC paper Fig. 2 and MMCL paper Fig. 1a)
%%-----------------------------------------------------------------
%
% In this example, we validate the MMCM algorithm using a homogeneous
% cubic domain. This simulation generates the results
% shown in Fig. 2 in the paper.
%
% The cubic domain has a dimension of 60x60x60 mm with optical properties
% mua=0.001, mus=1, n=1.0 and g=0.01. The analytical solution
% can be computed by the cwdiffusion (for CW solution) and
% tddiffusion (for time-domain solutions) functions in the
% package of Monte Carlo eXtreme (MCX) under the mcx/utils directory.
%
%%-----------------------------------------------------------------

%addpath('/Please/add/path/to/mcx/utils/')
addpath('../../matlab/');

clear cfg
cfg.nphoton=1e6;
cfg.seed=1648335518;
[cfg.node,cfg.elem]=genT5mesh(0:2:60,0:2:60,0:2:60);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1,30.2,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.prop=[0 0 1 1;0.005 1.0 0.01 1.0];
cfg.debuglevel='TP';
cfg.isreflect=0;
cfg.method='elem';

cfg.gpuid=1;

b1gpu=mmclab(cfg);
b1gpu=sum(b1gpu.data,2);

%%
cfg.gpuid=-1;

b1cpu=mmclab(cfg);
b1cpu=sum(b1cpu.data,2);


%%-----------------------------------------------------------------
%% add paths to the necessary toolboxes
%%-----------------------------------------------------------------

c0=299792458000;

twin=[cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend];
gates=length(twin);

[xi,yi]=meshgrid(0.5:1:59.5,0.5:1:59.5);

[cutpos1,cutvalue1,facedata1]=qmeshcut(cfg.elem(:,1:4),cfg.node,b1cpu,[0 30.5 0; 0 30.5 1; 1 30.5 0]);
[cutpos2,cutvalue2,facedata2]=qmeshcut(cfg.elem(:,1:4),cfg.node,b1gpu,[0 30.5 0; 0 30.5 1; 1 30.5 0]);
Phicpu=griddata(cutpos1(:,1),cutpos1(:,3),cutvalue1,xi,yi);
Phigpu=griddata(cutpos2(:,1),cutpos2(:,3),cutvalue2,xi,yi);


%% DMMC

[cfg.node,cfg.elem]=meshgrid6(0:60:60,0:60:60,0:60:60);
cfg.elem(:,1:4)=meshreorient(cfg.node,cfg.elem(:,1:4));

cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.method='grid';

cfg.gpuid=1;

flux=mmclab(cfg);

Phidgpu=sum(flux.data,4);

%% DMMC 
cfg.gpuid=-1;

flux=mmclab(cfg);
Phidcpu=sum(flux.data,4);

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.5
%%-----------------------------------------------------------------
figure
clines = 0:-0.5:-8;
[xi,yi]=meshgrid(0.5:1:59.5,0.5:1:59.5);
srcs=[30.1,30.2,0];
dets=[xi(:) 30.2*ones(size(xi(:))) yi(:)];

hold on
[c h2]=contourf(xi,yi, log10(max(squeeze(Phicpu*cfg.tstep),1e-8)), clines, 'k-', 'linewidth', 2);
contour(xi,yi,log10(abs(Phigpu*cfg.tstep)),clines,'g--','linewidth',2)
contour(xi,yi,log10(abs(squeeze(Phidcpu(1:end-1,31,1:end-1))'*cfg.tstep)),clines,'b--','linewidth',2)
contour(xi,yi,log10(abs(squeeze(Phidgpu(1:end-1,31,1:end-1))'*cfg.tstep)),clines,'r-','linewidth',1)

axis equal  
set(gca,'xlim',[0.5 59.5])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('B1-MMC','B1-MMCL','B1D-MMC','B1D-MMCL')
legend boxoff;
box on;
