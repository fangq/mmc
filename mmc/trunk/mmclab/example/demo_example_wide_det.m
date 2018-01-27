
addpath('../../matlab/');
addpath('../');
addpath('/home/mry09/Desktop/iso2mesh');

% wide-field illumination and 3 wide-field detection patterns for replay

clear cfg newcfg

cfg.nphoton=3e2;
cfg.seed=-1;
[cfg.node,cfg.face,cfg.elem]=meshabox([0 0 0], [60 60 20], 2, 2);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srctype='planar';
cfg.srcpos=[10,10,-1];
cfg.srcdir=[0 0 1];
cfg.srcparam1=[40 0 0 0];
cfg.srcparam2=[0 40 0 0];

cfg.elem = [cfg.elem, cfg.elemprop];

srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
    'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);

[cfg.node,cfg.elem] = mmcaddsrc(cfg.node,cfg.elem,srcdef);

detdef=struct('srctype','planar','srcpos',[10 10 20.1],'srcdir',[0 0 -1],...
    'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);

[cfg.node,cfg.elem] = mmcadddet(cfg.node, cfg.elem,detdef);

cfg.elemprop = cfg.elem(:,5);
cfg.elem = cfg.elem(:,1:4);

cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.prop=[0 0 1 1;0.005 1.0 0.01 1.0];
cfg.debuglevel='TP';
cfg.isreflect=0;
cfg.issaveexit=1;  % save detected photon exit position and angles
cfg.issaveseed=1;  % save detected photon seeds to replay later
% cfg.method=0;

newcfg=mmclab(cfg,'prep');  % preprocessing of the mesh to get the missing fields

% initial MC simulation
[cube, detp, ncfg, seeds]=mmclab(newcfg); % initial simulation

% prepare for replay (3 wide-field detection patterns simultaneously)

newcfg2 = newcfg;
newcfg2.outputtype = 'wl';
newcfg.replaydet = 0;
newcfg2.seed = seeds.data;
newcfg2.replayweight = mmcdetweight(detp.data,newcfg.prop);
newcfg2.replaytime = mmcdettime(detp.data,newcfg.prop);

newcfg2.detpos=[10. 10. 20.1 0
   10. 10. 20.1 0
   10. 10. 20.1 0];
newcfg2.detparam1=[40 0 0 2];
newcfg2.detparam2=[0 40 0 2];
pat=zeros(2,2,3);
pat(:,:,1)=[1,1;1,1];   % full detection pattern
pat(:,:,2)=[1,0;1,0];   % left half detection
pat(:,:,3)=[1,1;0,0];   % top half detection
newcfg2.detpattern=pat;
newcfg2.replaydetidx=mmcdetidx(detp.p,newcfg2.detpos,newcfg2.detparam1,newcfg2.detparam2);

% load newcfg2.mat;
[cube2,detp2,~,~]=mmclab(newcfg2);