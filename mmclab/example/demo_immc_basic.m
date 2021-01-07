% In this example, we show how to run implicit MMC (iMMC). The results are the
% same as three benchmarks in the paper: https://doi.org/10.1364/BOE.411898
%
% For iMMC, you need to edit cfg.elem and cfg.node. Extra data is needed to
% run iMMC. The data formats of cfg.elem and cfg.node for iMMC are shown 
% below:
%% edge-based iMMC (e-iMMC)
%                            elem | labeled edge | N/A | edge radius | N/A
% cfg.elem has 12 columns: [  1:4       5:6        7:8      9:10      11:12 ]
%                            node | node radius
% cfg.node has 4 columns:  [  1:3        4     ]
% 
% labeled edge: integer values from 0 to 6 indicate the 6 local edges in
%               the tetrahedral element. For example, 0->[1 2] represents
%               the 0th local edge connecting the 1st and 2nd column in
%               cfg.elem. Below is the mapping between local edge index
%               and its two connceted local nodes.
%               0->[1 2], 1->[1 3], 2->[1 4], 3->[2 3], 4->[2 4], 5->[3 4]
%               6->not labeled.
% edge radius:  7th and 8th columns are the radii of labeled edges in 5th
%               and 6th columns
% node radius:  the radii of the current node (used for node-based iMMC)
% N/A:          reserved for other purposes or future use

%% node-based iMMC (n-iMMC)
%                            elem | N/A
% cfg.elem has 12 columns: [  1:4   5:12]
%                            node | node radius
% cfg.node has 4 columns:  [  1:3        4     ]
% 
% node radius:  the radii of the current node

%% face-based iMMC (f-iMMC)
%                            elem | labeled face | face thickness
% cfg.elem has 12 columns: [ 1:4       5:8             9:12      ]
% 
% labeled face:    integer values from 0 to 3 indicate the 4 local faces
%                  represented by the local nodes. The mapping between faces
%                  and local node indices are:
%                  0->[4,1,2], 1->[4,2,3], 2->[3,1,4], 3->[2,1,3]
%                  For example, 0->[4,1,2] stands for the 0th face is
%                  represented by 4th, 1st, 2nd column in cfg.elem.
%                  The negative number indicates the element index current
%                  element is referring. Detail can be seen in the paper.
% face thickness:  the thickness of the faces in 5:8 columns

%% preparation
% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')

% 3. elem and node for iMMC
load immc

%% edge-based iMMC, benchmark B1
clear cfg

cfg.implicit = 1;    % turn on edge-based and node-based immc

cfg.nphoton=1e6;
cfg.node = node_eimmc;
cfg.elem = elem_eimmc;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[0.5 0.5 1];
cfg.srcdir=[0 0 -1];
cfg.prop=[0 0 1 1;0.0458 35.6541 0.9000 1.3700; 23.0543 9.3985 0.9000 1.3700];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.method = 'grid';
cfg.unitinmm = 0.01;
cfg.steps = [0.01 0.01 0.01];
cfg.isreflect = 1;
cfg.gpuid = -1;

% run edge-based iMMC
flux_eimmc=mmclab(cfg);

%% node-based iMMC, benchmark B2
cfg.node = node_nimmc;
cfg.elem = elem_nimmc;

% run node-based iMMC
flux_nimmc=mmclab(cfg);


%% face-based iMMC, benchmark B3
cfg.implicit = 2;
cfg.node = node_fimmc;
cfg.elem = elem_fimmc;

% run node-based iMMC
flux_fimmc=mmclab(cfg);

%% plot the results
figure,
subplot(131),imagesc(log10(rot90(squeeze(flux_eimmc.data(50,1:100,1:100)))))
subplot(132),imagesc(log10(rot90(squeeze(flux_nimmc.data(50,1:100,1:100)))))
subplot(133),imagesc(log10(rot90(squeeze(flux_fimmc.data(50,1:100,1:100)))))