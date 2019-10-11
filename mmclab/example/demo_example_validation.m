%%-----------------------------------------------------------------
%% validate MMC with a homogeneous cubic domain
%  (see MMC paper Fig. 2)
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

cfg.nphoton=3e6;
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

cube=mmclab(cfg);
%cube=mmclab(cfg,'sse'); % this is faster
cube=cube.data;

%%-----------------------------------------------------------------
%% add paths to the necessary toolboxes
%%-----------------------------------------------------------------

c0=299792458000;

twin=[cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend];
gates=length(twin);

cwcb=sum(cube,2);
[cutpos,cutvalue,facedata]=qmeshcut(cfg.elem(:,1:4),cfg.node,cwcb,[0 30.2 0; 0 30.2 1; 1 30.2 0]);

[xi,yi]=meshgrid(0:2:60,0:2:60);
vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

%%-----------------------------------------------------------------
%% plot the time-domain TPSF at [30 14 10]
%%-----------------------------------------------------------------

figure

srcpos=[30.1,30.2,0];
detpos=[30 14 10];

hold on
semilogy((1:gates)/10,tddiffusion(0.005, 1, c0, 0, srcpos, detpos,twin),'r');
idx=find(ismember(cfg.node,detpos,'rows'));
semilogy((1:gates)/10,cube(idx,:),'+');

set(gca,'fontsize',20)
xlabel('t (ns)')
ylabel('Fluence TPSF (1/mm^2)')
set(gca,'yscale','log')
legend('Diffusion','MMCM')
legend boxoff;
box on;

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure
clines = -1.5:-0.5:-8;
[xi,yi]=meshgrid(0:2:60,0:2:60);
srcs=[30.1,30.2,0];
dets=[xi(:) 30.2*ones(size(xi(:))) yi(:)];
phicw=reshape(cwdiffusion(0.005, 1.0, 0, srcs, dets),size(xi));

hold on
[c h2]=contour(xi,yi, log10(max(squeeze(phicw),1e-8)), -1.5:-0.5:-4, 'k-' );
contour(xi,yi,log10(abs(vi*cfg.tstep)),clines,'r:')

axis equal  
set(gca,'xlim',[1 60])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MMC')
legend boxoff;
box on;
