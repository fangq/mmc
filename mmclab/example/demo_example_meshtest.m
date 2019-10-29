%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Validate MMC in a heterogeneous domain: sphere in a cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this example, we validate MMCM algorithm using a sphere object
% in a homogeneous cubic background. This simulation generates the results
% shown in Fig. 3 in the paper.
% 
% The cubic domain has a dimension of 60x60x60 mm with optical properties
% mua=0.002, mus=1, n=1.37 and g=0.01. A sphere is centered at [30 30 30]mm
% with a radius of 10mm. The optical properties of the sphere is
% mua=0.05, mus=5, n=1.37 and g=0.9. The analytical solution, approximated
% by a sphere inside a infinite slab, can be computed by the
% sphdiffusionslab function in mmc/matlab/ directory (this has already
% been done, the result for plane y=30 is saved in the
% sphdiffsemiinf.mat file).
% 
% To validate MMCM, we generate an FE-mesh for the sphere and the cube.
% Three meshes are tested: mesh0: a coarse FE mesh with only 10,000 nodes,
% mesh1: a uniform dense FE mesh with 60000 nodes, and mesh2: a mesh
% with higher density around the sphere surface and near the source.
% The case mesh1 and mesh2 correspond to "MMCM Mesh 1" and "MMCM Mesh 2"
% in the paper (Table 1), respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation

% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')


% create a surface mesh for a 10mm radius sphere
[no,el]=meshasphere([30 30 30],10,1.0);

%%-----------------------------------------------------------------
%% create the common parameters for all 3 simulations
%%-----------------------------------------------------------------
clear cfg

cfg.nphoton=3e7;
cfg.seed=1648335518;
cfg.srcpos=[30.,30.,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.prop=[0 0 1 1;0.002 1.0 0.01 1.37;0.050 5.0 0.9 1.37];
cfg.debuglevel='TP';
cfg.isreflect=0;
cfg.method='elem';

cfg=[cfg cfg cfg];

%%-----------------------------------------------------------------
%% tetrahedral mesh generation
%%-----------------------------------------------------------------

% generate a coarse volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 20

ISO2MESH_SESSION='mmcmesh2_';

srcpos=[30. 30. 0.];
fixednodes=[30.,30.,0.05; 30 30 30];
nodesize=[ones(size(no,1),1) ; 0.5; 3];
nfull=[no;fixednodes];
[cfg(3).node,cfg(3).elem,face3]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,8,[30 30 30],[],[1.5 1.5 1.5 1.5 5 5 5 5]);
[cfg(3).node,cfg(3).elem]=sortmesh(srcpos,cfg(3).node,cfg(3).elem,1:4);
cfg(3).elem(:,5)=cfg(3).elem(:,5)+1;

% generate a dense volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 5

ISO2MESH_SESSION='mmcmesh1_';

nodesize=[1*ones(size(no,1),1) ; 1; 1];
[cfg(2).node,cfg(2).elem,face2]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,2,[30 30 30],[],[1 1 1 1 1 1 1 1]);
[cfg(2).node,cfg(2).elem]=sortmesh(srcpos,cfg(2).node,cfg(2).elem,1:4);
cfg(2).elem(:,5)=cfg(2).elem(:,5)+1;

% reduce the surface node numbers to 20%

ISO2MESH_SESSION='mmcmesh0_';

% using the coarse spherical surface, we generate a coarse volumetric
% mesh with maximum volume of 10

nodesize=[ones(size(no,1),1); 3];
nfull=[no; 30 30 30];
[cfg(1).node,cfg(1).elem,face1]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,10,[30 30 30],[],[2 2 2 2 5 5 5 5]);
[cfg(1).node,cfg(1).elem]=sortmesh(srcpos,cfg(1).node,cfg(1).elem,1:4);
cfg(1).elem(:,5)=cfg(1).elem(:,5)+1;

clear ISO2MESH_SESSION


%%-----------------------------------------------------------------
%% running MMC simulations for all 3 cases
%%-----------------------------------------------------------------

flux=mmclab(cfg);

%%-----------------------------------------------------------------
%% plotting the results
%%-----------------------------------------------------------------

c0=299792458000;

twin=[cfg(1).tstart+cfg(1).tstep/2:cfg(1).tstep:cfg(1).tend];
gates=length(twin);
clines=[-1:0.5:8]-10;

[xi,yi]=meshgrid(1:60,0:60);

%%-----------------------------------------------------------------
%% generate/load analytical solution for sphere inside infinite slab
%%-----------------------------------------------------------------

%[phi_ana,xa,ya,za]=sphdiffusionslab(0,0,60,-22:0.8:22,0,-30:0.8:10);
%save sphdiffsemiinf.mat phi_ana xa ya za

load ../../examples/meshtest/sphdiffsemiinf.mat
idx=find((xa(:)<-12 | xa(:)>12) & za(:)>-5);
phi_ana(idx)=nan;
idx=find((xa(:)<-10 | xa(:)>10) & za(:)>0);
phi_ana(idx)=nan;

%%-----------------------------------------------------------------
%% generate the contour of the inclusion
%%-----------------------------------------------------------------

[xcirc,ycirc] = cylinder(10,200);
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;

% create the voxel-contour of the sphere for MCX
%dim=60;
%[xv,yv,zv]=meshgrid(1:dim,1:dim,1:dim);
%dist=(xv-30).^2+(yv-30).^2+(zv-30).^2;

%vv=zeros(size(xv));
%vv(dist<100)=1;

%c=contourc(squeeze(vv(:,30,:)),1);
%plot(c(1,2:end),c(2,2:end),'c--')

%%-----------------------------------------------------------------
%% plot sphere 1
%%-----------------------------------------------------------------
node1=cfg(1).node;
elem1=cfg(1).elem;
mesh0=flux(1).data;
s1=sum(mesh0,2)*cfg(1).tstep;
[cutpos,cutvalue,facedata]=qmeshcut(elem1(:,1:4),node1,s1,[0 30 0; 0 30 1; 1 30 0]);
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])

vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana)),clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MMCM')
legend boxoff;
box on;

%%-----------------------------------------------------------------
%% plot sphere 2
%%-----------------------------------------------------------------

node2=cfg(2).node;
elem2=cfg(2).elem;
mesh1=flux(2).data;
s2=sum(mesh1,2)*cfg(1).tstep;
[cutpos,cutvalue,facedata]=qmeshcut(elem2(:,1:4),node2,s2,[0 30 0; 0 30 1; 1 30 0]);
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])
vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana)),clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MMCM')
legend boxoff;
box on;

%%-----------------------------------------------------------------
%% plot sphere 3
%%-----------------------------------------------------------------

node3=cfg(3).node;
elem3=cfg(3).elem;
mesh2=flux(3).data;
s3=sum(mesh2,2)*cfg(1).tstep;
[cutpos,cutvalue,facedata]=qmeshcut(elem3(:,1:4),node3,s3,[0 29 0; 0 29 1; 1 29 0]);

vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana)),clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MMCM')
legend boxoff;
box on;
