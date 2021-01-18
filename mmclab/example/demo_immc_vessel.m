% This example shows how to run e-iMMC on a vessel network

%% load vessel data
% node:   nodes of all vessel segments
% noder:  the radius corresponding to all nodes
% edge:   vessel segment index
% nbox:   nodes of bounding box
% fbox:   faces of bounding box

addpath(genpath('/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/'))
addpath('../../matlab')

load vessel

%% 1. Mesh generation

% (a) add the bounding box (nbox and fbox) into node and face
offset = size(nbox,1);
edge = edge+offset;  % change the node index after inserting bounding box
noder = [zeros(offset,1); noder];  % add offset into noder
node = [nbox; node];  % add nodes of bounding box into node
fedge = [edge(:,1) edge];  %  convert edge to face for mesh generation
face = [fbox; fedge];  %  add faces of bounding box into face

% (b) generate .poly file
generatePoly('vessel_mesh',node,face);

% (c) generate mesh files using tetgen
tic
% !/path/to/your/iso2mesh/bin/tetgen1.5.mexa64 -YY vessel_mesh.poly
!/space/neza/2/users/yaoshen/NEU/Research/iso2mesh/bin/tetgen1.5.mexa64 -YY vessel_mesh.poly
t=toc

% (d) read and plot vessel mesh
[nodeVessel, elemVessel, faceVessel]=readMesh('vessel_mesh');

%% 2. Label local edges

% 1. label the local edges which are vessels for all elements
[localEdgeIndex, edgeSet] = labelElements(edge,elemVessel);

% 2. split the elements wich over 2 edges labeled as vessel
indexForSplit = find(sum(logical(localEdgeIndex),2)>2);
if ~isempty(indexForSplit)
    [nodeVessel,elemVessel,localEdgeIndex] = relabelElement(elemVessel,localEdgeIndex,indexForSplit,nodeVessel,edgeSet);
    noder = [noder; zeros(size(nodeVessel,1)-size(noder,1),1)];
end

% 3. get radius for all edges (comment this line if you have your own vessel radii)
radii = getVesselRadii(elemVessel,localEdgeIndex,noder);

%% 3. Run mmc

clear cfg

cfg = inputMesh2Cfg(elemVessel,localEdgeIndex,radii,nodeVessel,noder);

res = 0.001;  % resolution: 0.001mm
cfg.steps = [0.001 0.001 0.001];
cfg.elem(:,9:12) = cfg.elem(:,9:12)*res;
cfg.node = cfg.node*res;
cfg.srctype = 'pencil';
cfg.nphoton=1e6;
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[150.5 150.5 605]*res;
cfg.srcdir=[0 0 -1];
cfg.prop=[0 0 1 1;
    0.0458   35.6541    0.9000    1.3700;   % dermis
    23.0543    9.3985    0.9000    1.3700]; % blood
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.method = 'grid';
cfg.gpuid = -1;
cfg.isreflect = 1;
cfg.outputtype = 'energy';

tic
flux=mmclab(cfg);
runtime = toc

figure,imagesc(log10(rot90(squeeze(flux.data(:,150,:)))))

%% %%%%%%%%%%% functions %%%%%%%%%%%

%% function: plot edges
function plotEdges(edge,node)
for i=1:size(edge,1)
    plot3([node(edge(i,1),1), node(edge(i,2),1)],...
        [node(edge(i,1),2), node(edge(i,2),2)],...
        [node(edge(i,1),3), node(edge(i,2),3)],'r');
    hold on
end
end

%% function: compute element centroid
function centroid = elemcentroid(node,elem)
% compute the centroid of all elements
centroid = zeros(size(elem,1),3);
for i=1:size(elem,1)
    n1 = node(elem(i,1),:);
    n2 = node(elem(i,2),:);
    n3 = node(elem(i,3),:);
    n4 = node(elem(i,4),:);
    n = [n1; n2; n3; n4];
    centroid(i,:) = mean(n);
end

end

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

%% function: label local edges for all elements
function [localEdgeIndex,edgeSet] = labelElements(edge,elem)
% label the local edges for all input element
n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
edge = sort(edge,2);          % sort to unify the edge representation
edgeSet = java.util.HashSet;  % set storing all vessel edges
for i=1:size(edge,1)
    edgeSet.add(num2str(edge(i,:)));
end

localEdgeIndex = zeros(size(elem,1),size(n2e,1));
fprintf('Getting local edge indices for each element...\n')
for i=1:size(elem,1)
   c = 1;
   for j=1:size(n2e,1)
       % check is local edge is in the edgeSet
       curLocalEdge = sort([elem(i,n2e(j,1)),elem(i,n2e(j,2))],2);
       if edgeSet.contains(num2str(curLocalEdge))
           localEdgeIndex(i,c) = j;
           c = c+1;
       end
   end
end
fprintf('Finished index extraction.\n')

end

%% function: relabel the elements
function [node,elemInserted,edgeIndexInserted] = relabelElement(elem,edgeIndex,indexForSplit,node,edgeSet)

fprintf('Relabeling elements with labeled edges larger than 2...\n')
faceorder = [1 2 3; 1 2 4; 1 3 4; 2 3 4];

% get elements that need to be split
elemSplit = elem(indexForSplit,:);

% get centroids and add them into nodes
centroid = elemcentroid(node,elemSplit);
offset = size(node,1);
node = [node; centroid];

%% get new elements and new labelings

% 4 new elements will be generated for each element to be split
splitElem = cell(length(indexForSplit),1);  % store new elements after spliting
splitLocalEdge = cell(length(indexForSplit),1);  % store local edge labeling for the new elements
for i=1:size(elemSplit,1)
    % one element can be split into 4 new elements
    newElem = zeros(size(faceorder,1),4);  % 4 nodes in an element
    newLocalEdge = zeros(size(faceorder,1),6);  % 6 edges in an element
    for j=1:size(faceorder,1)
        % offset+i is the node index of the centroid
        newElem(j,:) = [offset+i,...
                        elemSplit(i,faceorder(j,1)),...
                        elemSplit(i,faceorder(j,2)),...
                        elemSplit(i,faceorder(j,3))];
        newLocalEdge(j,:) = labelSingleElem(newElem(j,:),edgeSet);
        if sum(logical(newLocalEdge(j,:)))>2
           warning('Split failed')
        end
    end
    splitElem{i} = newElem;
    splitLocalEdge{i} = newLocalEdge;
end

%% Insert new elements to original elements

if indexForSplit(1)==1
    elemInserted = [];
    edgeIndexInserted = [];
else
    elemInserted = elem(1:indexForSplit(1)-1,:);
    edgeIndexInserted = edgeIndex(1:indexForSplit(1)-1,:);
end
elemInserted = [elemInserted; splitElem{1}];
edgeIndexInserted = [edgeIndexInserted; splitLocalEdge{1}];

for i=2:length(indexForSplit)
    elemInserted = [elemInserted; elem(indexForSplit(i-1)+1:indexForSplit(i)-1,:)];
    elemInserted = [elemInserted; splitElem{i}];
    
    edgeIndexInserted = [edgeIndexInserted; edgeIndex(indexForSplit(i-1)+1:indexForSplit(i)-1,:)];
    edgeIndexInserted = [edgeIndexInserted; splitLocalEdge{i}];
end

if indexForSplit(end)<size(elem,1)
    elemInserted = [elemInserted; elem(indexForSplit(end)+1:end,:)];
    edgeIndexInserted = [edgeIndexInserted; edgeIndex(indexForSplit(end)+1:end,:)];
end
fprintf('Finished relabeling.\n')

end

%% function: label local edges for single element
function localEdge = labelSingleElem(elem,edgeSet)
% label a single element
n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
localEdge = zeros(1,size(n2e,1));
c = 1;
for i=1:size(n2e,1)
    % check if the edge is in the set
    curEdge = sort([elem(n2e(i,1)) elem(n2e(i,2))],2);
    if edgeSet.contains(num2str(curEdge))
        localEdge(c) = i;
        c = c+1;
    end
end

end

%% function: obtain the radii for the labeled edges
function radii = getVesselRadii(elem,index,noder)
n2e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
radii = zeros(size(index,1),2);
for i=1:size(index,1)
    % edge 1
    if index(i,1)~=0
        ei1 = n2e(index(i,1),:);
        radii(i,1) = max([noder(elem(i,ei1(1))),noder(elem(i,ei1(2)))]);
    end
    
    % edge 2
    if index(i,2)~=0
        ei2 = n2e(index(i,2),:);
        radii(i,2) = max([noder(elem(i,ei2(1))),noder(elem(i,ei2(2)))]);
    end
end

end

%% function: build cfg for immc
function cfg = inputMesh2Cfg(elem,index,radii,node,noder)

% turn on iMMC
cfg.implicit = 1;

% post-processing to fit to iMMC input format
index = index(:,1:4);
index(index==0) = 7;
% 0 to 5 represet local edge labeled as vessel, 6 represent not labeled as
% vessel
index = index-1;

node = [node noder];
elem = [elem index radii zeros(size(radii,1),2)];

cfg.node = node;
cfg.elem = elem;
end