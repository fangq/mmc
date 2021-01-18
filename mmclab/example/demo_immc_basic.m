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

%% add path
% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');
addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')

%% edge-based iMMC, benchmark B1
% load immc

% (a) generate bounding box and insert edge
[nbox,ebox] = meshgrid6(0:1,0:1,0:1);
fbox = volface(ebox);
EPS = 0.001;
nbox = [nbox; [1-EPS 0.5 0.5]; [EPS 0.5 0.5]];  % insert new nodes (node 9 and 10)
fbox = [fbox; [9 9 10]];  % insert new edge coneected by node 9 and 10

% (b) generate .poly file
generatePoly('eimmc_mesh',nbox,fbox);

% (c) generate mesh files using tetgen
tic
% !/path/to/your/iso2mesh/bin/tetgen1.5.mexa64 -YY vessel_mesh.poly
!/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/bin/tetgen1.5.mexa64 -YY eimmc_mesh.poly
t=toc

% (d) read and plot vessel mesh
[node_eimmc, elem_eimmc, face_eimmc]=readMesh('eimmc_mesh');

% (e) label the edge that has node 9 and 10 and add radii
elem_eimmc = [elem_eimmc 6*ones(size(elem_eimmc))];
elem_eimmc([10 11 13 16],5) = [5 4 4 3];  % local edge index
elem_eimmc = [elem_eimmc zeros(size(elem_eimmc,1),4)];
elem_eimmc([10 11 13 16],9) = [0.1 0.1 0.1 0.1];  % edge radius = 0.1
node_eimmc = [node_eimmc zeros(size(node_eimmc,1),1)];
node_eimmc(9,4) = 0.1;  % add radius for inserted node

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
cfg.steps = [0.01 0.01 0.01];
cfg.isreflect = 1;
cfg.gpuid = -1;

% run edge-based iMMC
flux_eimmc=mmclab(cfg);

%% node-based iMMC, benchmark B2

% (a) generate bounding box and insert edge
[nbox,ebox] = meshgrid6(0:1,0:1,0:1);
fbox = volface(ebox);
nbox = [nbox; [0.5 0.5 0.5]];  % insert new nodes (node 9)

% (b) generate .poly file
generatePoly('nimmc_mesh',nbox,fbox);

% (c) generate mesh files using tetgen
tic
% !/path/to/your/iso2mesh/bin/tetgen1.5.mexa64 -YY vessel_mesh.poly
!/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/bin/tetgen1.5.mexa64 -YY nimmc_mesh.poly
t=toc

% (d) read and plot vessel mesh
[node_nimmc, elem_nimmc, face_nimmc]=readMesh('nimmc_mesh');

% (e) label the edge that has node 9 and 10 and add radii
elem_nimmc = [elem_nimmc 6*ones(size(elem_nimmc)) zeros(size(elem_nimmc))];
node_nimmc = [node_nimmc zeros(size(node_nimmc,1),1)];
node_nimmc(9,4) = 0.1;  % add radius for inserted node

% run node-based iMMC
cfg.elem = elem_nimmc;
cfg.node = node_nimmc;
flux_nimmc=mmclab(cfg);


%% face-based iMMC, benchmark B3

cfg.implicit = 2;  % set cfg.implicit to 2 to turn on f-iMMC

% (a) generate bounding box and insert edge
[nbox,ebox] = meshgrid6(0:1,0:1,0:1);
fbox = volface(ebox);

elem_fimmc = [ebox 6*ones(size(ebox)) zeros(size(ebox))];
elem_fimmc([3 4 5 6],5) = [-4 1 -6 1];
elem_fimmc([3 4 5 6],9) = [0.1 0.1 0.1 0.1];
node_fimmc = [nbox zeros(size(nbox,1),1)];

cfg.node = node_fimmc;
cfg.elem = elem_fimmc;

% run node-based iMMC
flux_fimmc=mmclab(cfg);

%% plot the results
figure,
subplot(131),imagesc(log10(rot90(squeeze(flux_eimmc.data(50,1:100,1:100)))))
subplot(132),imagesc(log10(rot90(squeeze(flux_nimmc.data(50,1:100,1:100)))))
subplot(133),imagesc(log10(rot90(squeeze(flux_fimmc.data(50,1:100,1:100)))))

%% function: generate .poly
function generatePoly(sessionid,node,face)
% generate the .poly file in current folder
fid=fopen([sessionid,'.poly'],'wt');
% === node ===
[nn,nd] = size(node);
fprintf(fid,'#node list\n');
fprintf(fid,'%d %d %d %d\n',nn,nd,0,0);
fprintf(fid,'%d %6.16f %6.16f %6.16f\n',[0:nn-1; node(:,1)'; node(:,2)'; node(:,3)']);

% === face ===
[ne,~] = size(face);
face = face-1;
fprintf(fid,'#face list\n');
fprintf(fid,'%d %d\n',ne,1);
for i=1:ne
    fprintf(fid,'%d %d\n',1,0);
    fprintf(fid,'%d %d %d %d\n',3,face(i,1),face(i,2),face(i,3));
end

% === hole ===
fprintf(fid,'#hole list\n');
fprintf(fid,'%d\n',0);

fclose(fid);
end

%% function: read generated .elem, .face and .node files
function [node_tet,elem_tet,face_tet]=readMesh(sessionid)
% read the .elem, .face and .node files generated from tetgen
% and write into elem_tet, face_tet and node_tet in MATLAB
fileID = fopen([sessionid '.1.ele'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
A = [A(1:3); 0; 0; A(4:end)];
elem_tet = reshape(A,[5 A(1)+1])';

fileID = fopen([sessionid '.1.face'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
A = [A(1:2); 0; 0; 0; A(3:end)];
face_tet = reshape(A,[5 A(1)+1])';

fileID = fopen([sessionid '.1.node'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
node_tet = reshape(A,[4 A(1)+1])';

elem_tet = elem_tet(2:end,2:end);
elem_tet = elem_tet+1;
face_tet = face_tet(2:end,2:end);
face_tet = face_tet+1;
node_tet = node_tet(2:end,2:end);

end