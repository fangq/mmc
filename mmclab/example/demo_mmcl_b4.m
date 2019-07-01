%% create the skin-vessel benchmark mesh
[no,fc]=latticegrid([0 200],[0 200],[20 32 200]); % create a 3-layer tissue
no(end,:)=no(end,:)+1e-5;

fc2=cell2mat(fc);
fc=[fc2(:,[1 2 3]); fc2(:,[1 3 4])];

[ncy,fcy]=meshacylinder([-1,99.5,99.5],[201,99.5,99.5],20,5); % add the vessel
[newnode,newelem]=surfboolean(no,fc,'first',ncy,fcy);  % merge the two domains

c0=[10,100,150,26]';
seeds=[ones(4,2)*100, c0];  % define the regions by index

[cfg.node,cfg.elem]=s2m(newnode,newelem(:,1:3),1,30,'tetgen',seeds,[],'-Y -A'); % creating the merged mesh domain

voxellen=0.005;
cfg.node=cfg.node*voxellen;
cfg.method='elem';

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
cfg.gpuid=1;
b4gpu=mmclab(cfg);
b4gpu=b4gpu.data;

%% saving coarse mesh
cfg.gpuid=-1;
b4cpu=mmclab(cfg);
b4cpu=b4cpu.data;

%% regenerate the mesh using fine mesh
[cfg.node,cfg.elem]=s2m(newnode,newelem(:,1:3),1,30,'tetgen',seeds,[]); % creating the merged mesh domain
cfg.node=cfg.node*voxellen;
cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);

cfg.method='grid';

[cfg.node,cfg.elem] = mmcaddsrc(cfg.node,[cfg.elem cfg.elemprop],...
    mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]));

%% saving coarse mesh
cfg.gpuid=1;
b4dgpu=mmclab(cfg);
b4dgpu=sum(b4dgpu,4);

%% saving coarse mesh
cfg.gpuid=-1;
b4dcpu=mmclab(cfg);
b4dcpu=sum(b4dcpu,4);

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure
clines = 0:-0.5:-8;
[xi,yi]=meshgrid(0:2:60,0:2:60);
srcs=[30.1,30.2,0];
dets=[xi(:) 30.2*ones(size(xi(:))) yi(:)];

hold on
[c h2]=contourf(xi,yi, log10(max(squeeze(Phicpu*cfg.tstep),1e-8)), clines, 'k-' );
contour(xi,yi,log10(abs(Phigpu*cfg.tstep)),clines,'r:')
contour(xi,yi,log10(abs(squeeze(b4dcpu(1:2:end,30,1:2:end))'*cfg.tstep)),clines,'b--')
contour(xi,yi,log10(abs(squeeze(b4dgpu(1:2:end,30,1:2:end))'*cfg.tstep)),clines,'y--')

axis equal  
set(gca,'xlim',[1 60])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('MMC','MMCL','D-MMC','D-MMCL')
legend boxoff;
box on;
