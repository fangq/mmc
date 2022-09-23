clear cfg

%% create surface mesh
[no_box,el_box]=meshgrid6(0:60:60,0:60:60,0:60:60);
fc_box=volface(el_box);

[no_sphere,fc_sphere]=meshasphere([30 30 30],25,3);
[no,fc]=mergemesh(no_box,fc_box,no_sphere,fc_sphere);

%% create volume mesh
ISO2MESH_TETGENOPT='-Y -A';
[cfg.node,cfg.elem]=surf2mesh(no,fc,[0 0 0],[60 60 60],1.0,100,...
    [1,1,1;30,30,30]);

%% set up cfg
cfg.nphoton=1e8;

cfg.srcpos=[30 30 0.01];
cfg.srcdir=[0 0 1];

cfg.prop=[0.000,  0, 1, 1;
          0.005,  1, 0, 1.37;  %box
          0.010, 10, 0, 1]; %sphere

% time-gate
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;

% dual grid MC
cfg.method='grid';

% output energy deposition
cfg.outputtype='fluence';

% gpu setting
cfg.gpuid=1;

% save configuration
mmc2json(cfg,'optix');