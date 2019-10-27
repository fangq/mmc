%% create the skin-vessel benchmark mesh
[no,fc]=latticegrid([0 200],[0 200],[20 32 200]); % create a 3-layer tissue
no(end,:)=no(end,:)+1e-5;

fc2=cell2mat(fc);
fc=[fc2(:,[1 2 3]); fc2(:,[1 3 4])];

[ncy,fcy]=meshacylinder([-1,99.5,99.5],[201,99.5,99.5],20,5); % add the vessel
[newnode,newelem]=surfboolean(no,fc,'first',ncy,fcy);  % merge the two domains

c0=[10,100,150,26]';
seeds=[ones(4,2)*100, c0];  % define the regions by index

[cfg.node,cfg.elem]=s2m(newnode,newelem(:,1:3),1,30,'tetgen',seeds,[],'-YY -A'); % creating the merged mesh domain

voxellen=0.005;
cfg.node=cfg.node*voxellen;
cfg.method='g';

figure; 
subplot(121);
plotmesh(cfg.node,cfg.elem);

cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);

%% define other properties

cfg.prop=[0.0000         0.0    1.0000    1
    3.5640e-05    1.0000    1.0000    1.3700
   23.0543    9.3985    0.9000    1.3700
    0.0458   35.6541    0.9000    1.3700
    1.6572   37.5940    0.9000    1.3700];

cfg.srcpos=[100 100 -1]*voxellen;
cfg.srcdir=[0 0 1];

cfg.tstart=0;
cfg.tend=5e-8;
cfg.tstep=5e-8;
% cfg.outputtype='energy'; %energy deposition in mmc varys with elem volume
cfg.outputtype='flux';
cfg.minenergy=0.01;

cfg.srctype='disk';
cfg.srcparam1=[0.3 0 0 0];

%% define wide-field disk source by extending the mesh to the widefield src
srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
    'srcparam1',cfg.srcparam1,'srcparam2',[]);

[cfg.node,cfg.elem] = mmcaddsrc(cfg.node,[cfg.elem cfg.elemprop],...
    mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]));


cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);

%% other simulation information

cfg.nphoton=1e7;
cfg.seed=1648335518;

cfg.debuglevel='TP';
cfg.isreflect=0;

%% saving coarse mesh
mmc2json(cfg,'dmmc_skinvessel.json');

%% regenerate the mesh using fine mesh
[cfg.node,cfg.elem]=s2m(newnode,newelem(:,1:3),1,30,'tetgen',seeds,[]); % creating the merged mesh domain
cfg.node=cfg.node*voxellen;
cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);

[cfg.node,cfg.elem] = mmcaddsrc(cfg.node,[cfg.elem cfg.elemprop],...
    mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]));

mmc2json(cfg,'skinvessel.json');

