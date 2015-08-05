%%-----------------------------------------------------------------
%% Comparing widefield MMC and MCX in a homogeneous cubic domain
%  (see MMC paper Fig. 2)
%%-----------------------------------------------------------------
%
% In this example, we compare MMC and MCX using both simple and
% widefield sources.
%
% The cubic domain has a dimension of 60x60x60 mm with optical properties
% mua=0.001, mus=1, n=1.0 and g=0.01. The analytical solution for 
% pencilbeam source can be computed by the cwdiffusion (for CW solution) 
% and tddiffusion (for time-domain solutions) functions in the
% package of Monte Carlo eXtreme (MCX) under the mcx/utils directory.
%
%%-----------------------------------------------------------------

%addpath('/Please/add/path/to/mcx/utils/')
addpath('../../matlab/');
addpath('../../mmclab/');

%%-----------------------------------------------------------------
%% run mmclab for the 60x60x60 homogeneous cubic domain
%%-----------------------------------------------------------------

planarsrc=1; % set to 1 to compare planar src, set to 0 for pencil beam

clear cfg;

cfg.nphoton=3e7;
cfg.seed=27182818;
[cfg.node,tmp,cfg.elem]=meshabox([0 0 0],[60 60 60], 2*sqrt(2), 2^3/6);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1,30.2,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.prop=[0 0 1 1;0.005 1.0 0.01 1.0];
cfg.debuglevel='TP';
cfg.isreflect=0;

% define wide-field planar source

if(planarsrc)   % when planar src is used, re-tessellate the mesh
    cfg.srctype='planar';
    cfg.srcpos=[10 10 -1];
    cfg.srcparam1=[40 0 0 0];
    cfg.srcparam2=[0 40 0 0];

    % the following block before "end" can be removed as mmclab
    % has this built-in
    srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
                  'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);
    cfg.elem(:,5)=1;
    [cfg.node,cfg.elem] = mmcaddsrc(cfg.node,cfg.elem,...
         mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]));
    cfg.elemprop=cfg.elem(:,5);
    cfg.elem=cfg.elem(:,1:4);
end

%%-----------------------------------------------------------------
%% run mmclab
%%-----------------------------------------------------------------

cube=mmclab(cfg);
%cube=mmclab(cfg,'sse'); % this is faster
cube=cube.data;

%%-----------------------------------------------------------------
%% run mcxlab for the same domain
%%-----------------------------------------------------------------

cfgx=cfg;
cfgx=rmfield(cfgx,{'node','elem','elemprop','debuglevel'});
dim=60;
cfgx.vol=ones(dim,dim,dim);
cfgx.vol=uint8(cfgx.vol);
cfgx.nphoton=numel(cfgx.vol)/size(cube,1)*cfg.nphoton; % for mcx, nphoton 
                                                       % needs to be bigger 
                                                       % to match the noise 
                                                       % of mmc
cfgx.autopilot=1;
cfgx.gpuid=1;
cfgx.issrcfrom0=1; % treating the lower-outer corner of the grid as [0 0 0]
                   %to match mmc

cubex=mcxlab(cfgx);
cubex=cubex.data;

%%-----------------------------------------------------------------
%% Make a cross-section
%%-----------------------------------------------------------------

c0=299792458000;

twin=[cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend];
gates=length(twin);

cwcb=sum(cube,2);
cwcbx=sum(cubex,4);

[cutpos,cutvalue,facedata]=qmeshcut(cfg.elem(:,1:4),cfg.node,cwcb,'y=30.5');
                                         % slice at y=30.5 to offset half
                                         % grid to match mcx voxels

[xi,yi]=meshgrid(0.5:59.5,0.5:59.5);     % interpolate mmc solution at 
                                         % half-grid because mcx's voxel 
                                         % readings are at the center of 
                                         % the voxels
vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);
vix=squeeze(cwcbx(:,30,:))';

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.5
%%-----------------------------------------------------------------
figure
hold on
clines = -1.5:-0.5:-8;
if(~planarsrc)
    srcs=[30.1,30.2,0];
    dets=[xi(:) 30.2*ones(size(xi(:))) yi(:)];
    phicw=reshape(cwdiffusion(0.005, 1.0, 0, srcs, dets),size(xi));
    contour(xi,yi, log10(max(squeeze(phicw),1e-8)), -1.5:-0.5:-4, 'k-' ); 
                 % analytical is semi-infinite, only take the center part
end
contour(xi,yi,log10(abs(vi*cfg.tstep)),clines,'r:')
contour(xi,yi,log10(abs(vix*cfg.tstep)),clines,'b--')

axis equal
set(gca,'xlim',[1 60]);
set(gca,'fontsize',20);
xlabel('x (mm)');
ylabel('z (mm)');
if(planarsrc)
    legend('MMC','MCX');
else
    legend('Diffusion','MMC','MCX');
end
legend boxoff;
box on;
