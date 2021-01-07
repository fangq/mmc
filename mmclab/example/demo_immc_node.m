% Temporary path
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))

load immc

clear cfg

cfg.implicit = 1;  % turn on node-based immc

cfg.nphoton=1e6;
cfg.node = node_nimmc;
cfg.elem = elem_nimmc;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[0.5 0.5 1];
cfg.srcdir=[0 0 -1];
cfg.prop=[0 0 1 1;0.0458 35.6541 0.9000 1.3700;23.0543 9.3985 0.9000 1.3700];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.method = 'grid';
cfg.unitinmm = 0.01;
cfg.steps = [0.01 0.01 0.01];
cfg.isreflect=1;
cfg.gpuid = -1;

%% run node-based iMMC
flux=mmclab(cfg);

%% plot the result
figure,imagesc(log10(rot90(squeeze(flux.data(50,1:100,1:100)))))