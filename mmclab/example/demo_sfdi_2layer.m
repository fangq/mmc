%%-----------------------------------------------------------------
%% Simulating an SFDI source with a 2-layer brain model
%%-----------------------------------------------------------------
%
% In this example, we simulate an spatial-frequency domain imaging source
% using a 2-layer brain model.
%
% The cubic domain has a dimension of 60x60x30 mm with the 5mm top layer
% simulating the skull/scalp while the bottom 25mm layer simulating gray
% and white matters. The sample SFDI covers a 40 x 40 mm area from the top,
% and has kx=3 with a pi/3 phase offset in the x-direction.
%
%%-----------------------------------------------------------------

clear cfg;
clear all;

%% 2-layer mesh model

layercount=3;

if(layercount==2)
    [node,face,c0]=latticegrid([0 60],[0 60],[0 25 30]);
    c0(:,4)=[2;1];   % maximum element size for bottom (label 1) and top (label 2) layers
    
    % simulate the optical properties of skull and gray-matter of a human brain model
    % see http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
    cfg.prop=[0 0 1 1;0.02 9.0, 0.89 1.37;0.019 7.8 0.89 1.37]; 
else
    [node,face,c0]=latticegrid([0 60],[0 60],[0 20 25 30]); % if you like a 3-layer model
    c0(:,4)=[2;2;1];
    cfg.prop=[0 0 1 1;0.02 9.0, 0.89 1.37;0.004 0.009, 0.89 1.37;0.019 7.8 0.89 1.37]; 
end
[cfg.node,cfg.elem]=surf2mesh(node,face,[],[],1,[],c0);

figure; 
subplot(121);
plotmesh(cfg.node,cfg.elem);

cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);


%% add source and retessellate mesh

cfg.srctype='fourier';      % define an SFDI source
cfg.srcpos=[10 10 35];      % one corner of the illumination area
kx=3;                       % wave number in the x-dir
ky=0;                       % wave number in the x-dir
xphase=pi/3;                % phase offset in the x-dir, must < 2pi
yphase=0;                   % phase offset in the x-dir, must < 2pi
cfg.srcparam1=[40 0 0 kx+xphase/(2*pi)];   % kx is k-number in x direction
cfg.srcparam2=[0 40 0 ky+yphase/(2*pi)];
cfg.srcdir=[0 0 -1];

%% line 24-31 could possibly be deleted and replaced by built-in one-step command
srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
    'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);

[cfg.node,cfg.elem] = mmcaddsrc(cfg.node,[cfg.elem cfg.elemprop],...
    mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]));

cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);

%% other simulation information

cfg.nphoton=3e6;
cfg.seed=1648335518;

cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;

cfg.debuglevel='TP';
cfg.isreflect=1;
cfg.detpos=[30 30 0 2]; % detector position
cfg.method='elem';

%% mmc simulation

layer=mmclab(cfg);
%cube=mmclab(cfg,'sse'); % this is faster
layer=layer.data;
layercw=sum(layer,2);

% plot simulated photon profiles

subplot(122);
hold on;
qmeshcut(cfg.elem(cfg.elemprop>0,1:4),cfg.node,log10(layercw),'y=30','linestyle','none'); view(3)
qmeshcut(cfg.elem(cfg.elemprop>0,1:4),cfg.node,log10(layercw),'z=27','linestyle','none'); 
if(layercount~=2)
    qmeshcut(cfg.elem(cfg.elemprop>0,1:4),cfg.node,log10(layercw),'z=22','linestyle','none'); 
end
box on;
axis equal
view(-56, 22);
