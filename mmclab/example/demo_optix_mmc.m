clear
clc
close all

%% create surface mesh
[no_box,el_box]=meshgrid6(0:60:60,0:60:60,0:60:60);
fc_box=volface(el_box);

[no_sphere,fc_sphere]=meshasphere([30 30 30],25,3);
[no,fc]=mergemesh(no_box,fc_box,no_sphere,fc_sphere);

%% create volume mesh
ISO2MESH_TETGENOPT='-Y -A';
[node,elem,face]=surf2mesh(no,fc,[0 0 0],[60 60 60],1.0,100,...
    [1,1,1;30,30,30]);

%% set up cfg
clear cfg
cfg.nphoton=1e8;

cfg.node=node;
cfg.elem=elem;

cfg.srcpos=[30 30 0.01];
cfg.srcdir=[0 0 1];

cfg.prop=[0.000,  0, 1, 1;
          0.005,  1, 0, 1.37;  %box
          0.010, 10, 0, 1]; %sphere

% disable reflection
cfg.isreflect=1;

% time-gate
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;

% dual grid MC
cfg.method='grid';

% output energy deposition
cfg.outputtype='energy';

% disable normalization
cfg.isnormalized=0;

% gpu setting
cfg.gpuid=1;

%% run optix MMC simulation
cfg.compute='optix';

output=mmclab(cfg);
energy=output.data(1:end-1,1:end-1,1:end-1,:);
energyoptix=sum(energy,4);

%% run opencl MMC simulation
cfg.compute='opencl';
output=mmclab(cfg);
energy=output.data(1:end-1,1:end-1,1:end-1,:);
energyopencl=sum(energy,4);

%% compare results
figure;
clines=-20:2:10;
contourf(log(squeeze(energyopencl(30,:,:))'),clines,'k-','displayname','OpenCL');
hold on;
contour(log(squeeze(energyoptix(30,:,:))'),clines,'r--','displayname','Optix');
axis equal;
legend;
colorbar;
title("Energy deposition (10^8) photons");