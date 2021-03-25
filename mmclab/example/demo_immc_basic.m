% In this example, we show how to run implicit MMC (iMMC). The results are the
% same as three benchmarks in the paper: https://doi.org/10.1364/BOE.411898
%
% For iMMC, you need to edit cfg.elem and cfg.node. Extra data is needed to
% run iMMC. The data formats of cfg.elem and cfg.node for iMMC are shown 
% below:
%% edge-based iMMC (e-iMMC)
%
% cfg.edgeroi - a 6-column array, storing the cylindrial edge-roi radii in the nchoosek order
%               [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]. A 0 indicates no ROI. Up to 2 edges are supported
%
% cfg.noderoi - a single column/row vector, storing the spherical node-roi radii
%
%% node-based iMMC (n-iMMC)
%
% cfg.noderoi - a single column/row vector, storing the spherical node-roi radii
%
%% face-based iMMC (f-iMMC)
%
% cfg.faceroi - a 4 column array, storing the slab-thickness of face-roi in the nchoosek order
%               [(1,2,3),(1,2,4),(1,3,4),(2,3,4)]. A 0 indicates no ROI.

%% add path
% 1. you need to add the path to iso2mesh toolbox if not already
%addpath('/path/to/iso2mesh/toolbox/');

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')

%% edge-based iMMC, benchmark B1

% (a) generate bounding box and insert edge
[nbox,ebox]=meshgrid6(0:1,0:1,0:1);
fbox=volface(ebox);
EPS=0.001;
nbox=[nbox; [1-EPS 0.5 0.5]; [EPS 0.5 0.5]];  % insert new nodes (node 9 and 10)
fbox=[fbox; [9 9 10]];  % insert new edge coneected by node 9 and 10

clear cfg0 cfg;

cfg0.nphoton=1e6;
cfg0.srcpos=[0.5 0.5 1];
cfg0.srcdir=[0 0 -1];
cfg0.prop=[0 0 1 1;0.0458 35.6541 0.9000 1.3700; 23.0543 9.3985 0.9000 1.3700];
cfg0.tstart=0;
cfg0.tend=5e-9;
cfg0.tstep=5e-9;
cfg0.debuglevel='TP';
cfg0.method='grid';
cfg0.steps=[0.01 0.01 0.01];
cfg0.isreflect=1;
cfg0.gpuid=-1;


cfg=cfg0;

% (b) generate mesh
[cfg.node,cfg.elem]=s2m(nbox,num2cell(fbox,2),1,100,'tetgen1.5',[],[],'-YY');
cfg.elemprop=ones(size(cfg.elem,1),1);

% (c) label the edge that has node 9 and 10 and add radii

% find all elements that contains special edge between nodes 9-10
edgeid=reshape(ismember(sort(meshedge(cfg.elem),2),[9 10],'rows'),[size(cfg.elem,1),6]);

% create the edgeroi input - every appearance of that edge, we define its radius
cfg.edgeroi=zeros(size(cfg.elem,1),6);
cfg.edgeroi(edgeid)=0.1;

cfg.noderoi=zeros(size(cfg.node,1),1);
cfg.noderoi(9)=0.1;


% run edge-based iMMC
flux_eimmc=mmclab(cfg);

%% node-based iMMC, benchmark B2

% (a) generate bounding box and insert edge
[nbox,ebox]=meshgrid6(0:1,0:1,0:1);
fbox=volface(ebox);
nbox=[nbox; [0.5 0.5 0.5]];  % insert new nodes (node 9)

cfg=cfg0;

% (b) generate mesh
[cfg.node,cfg.elem]=s2m(nbox,num2cell(fbox,2),1,100,'tetgen1.5',[],[],'-YY');
cfg.elemprop=ones(size(cfg.elem,1),1);

% (c) label the edge that has node 9 and 10 and add radii
cfg.noderoi=zeros(size(cfg.node,1),1);
cfg.noderoi(9)=0.1;

% run node-based iMMC
flux_nimmc=mmclab(cfg);


%% face-based iMMC, benchmark B3

cfg=cfg0;

% (a) generate bounding box and insert edge
[cfg.node,cfg.elem]=meshgrid6(0:1,0:1,0:1);
cfg.elemprop=ones(size(cfg.elem,1),1);

% generate all triangular faces, label the top two triangles at z=1 as
% faceroi and define a h=0.1mm thin-membrane
facelist=nchoosek(1:4,3);
allfaces=[];
for i=1:length(facelist)
   allfaces=[allfaces;cfg.elem(:,facelist(i,:))];
end

% sum(allfaces>4,2)==3 gives the triangles that are made of nodes 5/6/7/8
cfg.faceroi=reshape(double(sum(allfaces>4,2)==3), size(cfg.elem,1),size(facelist,1));
cfg.faceroi(cfg.faceroi>0)=0.1;

% run node-based iMMC
flux_fimmc=mmclab(cfg);

%% plot the results
figure,
subplot(131)
imagesc(log10(rot90(squeeze(flux_eimmc.data(50,1:100,1:100)))))
title('edge-iMMC');axis equal
subplot(132)
imagesc(log10(rot90(squeeze(flux_nimmc.data(50,1:100,1:100)))))
title('node-iMMC');axis equal
subplot(133)
imagesc(log10(rot90(squeeze(abs(flux_fimmc.data(50,1:100,1:100))))))
title('face-iMMC');axis equal
