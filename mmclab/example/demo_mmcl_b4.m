%% create the skin-vessel benchmark mesh
[no,fc]=latticegrid([0 200],[0 200],[20 32 200]); % create a 3-layer tissue
no(end,:)=no(end,:);

fc2=cell2mat(fc);
fc=[fc2(:,[1 2 3]); fc2(:,[1 3 4])];

[ncy,fcy]=meshacylinder([-1,99.5,99.5],[201,99.5,99.5],20,5); % add the vessel
[newnode,newelem]=surfboolean(no,fc,'first',ncy,fcy);  % merge the two domains

c0=[10,100,150,26]';
seeds=[ones(4,2)*100, c0];  % define the regions by index

[cfg_mmc.node,cfg_mmc.elem]=s2m(newnode,newelem(:,1:3),1,4,'tetgen',seeds,[]); % creating the merged mesh domain

voxellen=0.005;
cfg_mmc.node=cfg_mmc.node*voxellen;
cfg_mmc.method='elem';

figure; 
subplot(121);
plotmesh(cfg_mmc.node,cfg_mmc.elem);

cfg_mmc.elemprop=cfg_mmc.elem(:,5);
cfg_mmc.elem=cfg_mmc.elem(:,1:4);

%% define other properties

cfg_mmc.prop=[0.0000         0.0    1.0000    1
    3.5640e-05    1.0000    1.0000    1.3700
   23.0543    9.3985    0.9000    1.3700
    0.0458   35.6541    0.9000    1.3700
    1.6572   37.5940    0.9000    1.3700];

cfg_mmc.srcpos=[100 100 -3]*voxellen;
cfg_mmc.srcdir=[0 0 1];

cfg_mmc.tstart=0;
cfg_mmc.tend=5e-8;
cfg_mmc.tstep=5e-8;
% cfg.outputtype='energy'; %energy deposition in mmc varys with elem volume
cfg_mmc.outputtype='flux';
cfg_mmc.minenergy=0.01;

cfg_mmc.srctype='disk';
cfg_mmc.srcparam1=[0.3 0 0 0];

%% define wide-field disk source by extending the mesh to the widefield src
srcdef=struct('srctype',cfg_mmc.srctype,'srcpos',cfg_mmc.srcpos,'srcdir',cfg_mmc.srcdir,...
    'srcparam1',cfg_mmc.srcparam1,'srcparam2',[]);

[cfg_mmc.node,cfg_mmc.elem] = mmcaddsrc(cfg_mmc.node,[cfg_mmc.elem cfg_mmc.elemprop],...
    mmcsrcdomain(srcdef,[min(cfg_mmc.node);max(cfg_mmc.node)]));


cfg_mmc.elemprop=cfg_mmc.elem(:,5);
cfg_mmc.elem=cfg_mmc.elem(:,1:4);

%% other simulation information

cfg_mmc.nphoton=1e8;
cfg_mmc.seed=1648335518;

cfg_mmc.debuglevel='TP';
cfg_mmc.isreflect=0;

%% saving coarse mesh
cfg_mmc.gpuid=1;
output1=mmclab(cfg_mmc);
b4gpu=output1.data*cfg_mmc.tstep;

%% saving coarse mesh
cfg_mmc.gpuid=-1;
output2=mmclab(cfg_mmc);
b4cpu=output2.data*cfg_mmc.tstep;

%% regenerate the mesh using fine mesh
% cfg_dmmc=cfg_mmc;
% [cfg_dmmc.node,cfg_dmmc.elem]=s2m(newnode,newelem(:,1:3),1,30,'tetgen',seeds,[],'-Y -A'); % creating the merged mesh domain
% cfg_dmmc.node=cfg_dmmc.node*voxellen;
% cfg_dmmc.elemprop=cfg_dmmc.elem(:,5);
% cfg_dmmc.elem=cfg_dmmc.elem(:,1:4);
% 
% cfg_dmmc.method='grid';
% cfg_dmmc.unitinmm=voxellen;
% 
% [cfg_dmmc.node,cfg_dmmc.elem] = mmcaddsrc(cfg_dmmc.node,[cfg_dmmc.elem cfg_dmmc.elemprop],...
%     mmcsrcdomain(srcdef,[min(cfg_dmmc.node);max(cfg_dmmc.node)]));
% 
% cfg_dmmc.elemprop=cfg_dmmc.elem(:,5);
% cfg_dmmc.elem=cfg_dmmc.elem(:,1:4);
% %% saving coarse mesh
% % cfg_dmmc.nphoton=9;
% % cfg_dmmc.gpuid=1;
% % b4dgpu=mmclab(cfg_dmmc);
% % b4dgpu=sum(b4dgpu.data,4);
% 
% %% saving coarse mesh
% cfg_dmmc.gpuid=-1;
% b4dcpu=mmclab(cfg_dmmc);
% b4dcpu=sum(b4dcpu.data,4);

%% interpolation
% x_min=min(cfg.node(:,1)); x_max=max(cfg.node(:,1)); bb_size(1,1)=floor(x_max-x_min);
% y_min=min(cfg.node(:,2)); y_max=max(cfg.node(:,2)); bb_size(1,2)=floor(y_max-y_min);
% z_min=min(cfg.node(:,3)); z_max=max(cfg.node(:,3)); bb_size(1,3)=floor(z_max-z_min);
[yi,zi]=ndgrid(0.0025:0.005:0.9975,0.1025:0.005:0.9975);
[cutpos1,cutvalue1,facedata1]=qmeshcut(cfg_mmc.elem(:,1:4),cfg_mmc.node,b4cpu,[0.4975 0 0;0.4975 0 1;0.4975 1 0]);
[cutpos2,cutvalue2,facedata2]=qmeshcut(cfg_mmc.elem(:,1:4),cfg_mmc.node,b4gpu,[0.4975 0 0;0.4975 0 1;0.4975 1 0]);
b4cpu_interpolated=griddata(cutpos1(:,2),cutpos1(:,3),cutvalue1,yi,zi);
b4gpu_interpolated=griddata(cutpos2(:,2),cutpos2(:,3),cutvalue2,yi,zi);

%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure
clines = -2:0.25:1;

hold on
[yyy,zzz]=ndgrid(-0.4975:0.005:0.4975,0.1025:0.005:0.9975);
contourf(yyy,zzz,log10(b4cpu_interpolated),clines,'k-','linewidth',2);colorbar
contour(yyy,zzz,log10(b4gpu_interpolated),clines,'m:','linecolor',[0.9100    0.4100    0.1700],'linewidth',2);

[xcirc,ycirc] = cylinder(0.1,10);
xcirc=xcirc(1,:)-0.0025;
ycirc=ycirc(1,:)+0.4975;
plot(xcirc,ycirc,'k--','linewidth',2);

x=[-0.5,0.5]; y=[0.16 0.16];
plot(x,y,'k--','linewidth',2);

%title('fluence \phi [W/mm^2/W.delivered]')
axis equal  
xlabel('y (mm)');xlim([-0.5 0.5]);
ylabel('z (mm)');ylim([0.1 1]);
axis equal
%legend('MMC','MMCL')
set(gca,'fontsize',20);
legend boxoff;
box on;

%% in case if you want to be consistant with Dr. Jacques' results, which is in unit: W/cm^2/W
figure
clines = 0:0.25:3;
hold on
[yyy,zzz]=ndgrid(-0.04975:0.0005:0.04975,0.01025:0.0005:0.09975);
contourf(yyy,zzz,log10(b4cpu_interpolated*100),clines,'k-','linewidth',2);colorbar
contour(yyy,zzz,log10(b4gpu_interpolated*100),clines,'g--','linewidth',2);

[xcirc,ycirc] = cylinder(0.01,10);
xcirc=xcirc(1,:)-0.00025;
ycirc=ycirc(1,:)+0.04975;
plot(xcirc,ycirc,'k--','linewidth',2);

x=[-0.05,0.05]; y=[0.016 0.016];
plot(x,y,'k--','linewidth',2);


title('fluence \phi [W/cm^2/W.delivered]')
xlabel('y [cm]');xlim([-0.05 0.05]);
ylabel('z [cm]');ylim([0.01 0.1]);
axis equal
legend('MMC','MMCL')
set(gca,'fontsize',20);
legend boxoff;
box on;