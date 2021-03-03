function cfg=demo_immc_vessel(cfg)
% This example shows how to run e-iMMC on a vessel network

%% load vessel data
% node:   nodes of all vessel segments
% noder:  the radius corresponding to all nodes
% edge:   vessel segment index
% nbox:   nodes of bounding box
% fbox:   faces of bounding box

%addpath('/path/to/iso2mesh/toolbox/');

if(nargin==0)

    vs=load('vessel.mat');

    edge=vs.edge;   % vessel edges
    fbox=vs.fbox;   % bounding box faces
    nbox=vs.nbox;   % bounding box nodes
    node=vs.node;   % nodes on the vessel network
    noder=vs.noder; % nodal radius of each node

    clear vs;

    %% 1. Mesh generation

    % (a) add the bounding box (nbox and fbox) into node and face
    offset=size(nbox,1);
    edge=edge+offset;  % change the node index after inserting bounding box
    noder=[zeros(offset,1); noder];  % add offset into noder
    node=[nbox; node];  % add nodes of bounding box into node
    fedge=[edge(:,1) edge];  %  convert edge to face for mesh generation
    face=[fbox; fedge];  %  add faces of bounding box into face

    % (b) generate .poly file
    [nodeVessel, elemVessel]=s2m(node,num2cell(face,2),1,100,'tetgen1.5',[],[],'-YY');

    %% 2. Label local edges

    % 1. label the local edges which are vessels for all elements
    [localEdgeIndex, edgeSet]=labelElements(edge,elemVessel);

    % 2. split the elements wich over 2 edges labeled as vessel
    indexForSplit=find(sum(logical(localEdgeIndex),2)>2);
    if ~isempty(indexForSplit)
        [nodeVessel,elemVessel,localEdgeIndex]=relabelElement(elemVessel,localEdgeIndex,indexForSplit,nodeVessel,edgeSet);
        noder=[noder; zeros(size(nodeVessel,1)-size(noder,1),1)];
    end

    % 3. get radius for all edges (comment this line if you have your own vessel radii)
    radii=getVesselRadii(elemVessel,localEdgeIndex,noder);

    %% 3. Run mmc

    cfg=inputMesh2Cfg(elemVessel,localEdgeIndex,radii,nodeVessel,noder);

    res=0.001;  % resolution: 0.001mm
    cfg.steps=[0.001 0.001 0.001];
    cfg.elem(:,9:12)=cfg.elem(:,9:12)*res;
    cfg.node=cfg.node*res;
    cfg.srctype='pencil';
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
    cfg.method='grid';
    cfg.gpuid=-1;
    cfg.isreflect=1;
    cfg.outputtype='energy';

end

tic
flux=mmclab(cfg);
toc

%% plotting
figure,imagesc(log10(rot90(squeeze(flux.data(:,150,:)))))
axis equal;
title('energy deposition map');


%% %%%%%%%%%%% functions %%%%%%%%%%%

%% function: label local edges for all elements
function [localEdgeIndex,edgeSet]=labelElements(edge,elem)
% label the local edges for all input element
n2e=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
edge=sort(edge,2);          % sort to unify the edge representation
edgeSet=java.util.HashSet;  % set storing all vessel edges
for i=1:size(edge,1)
    edgeSet.add(num2str(edge(i,:)));
end

localEdgeIndex=zeros(size(elem,1),size(n2e,1));
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


%% function: relabel the elements
function [node,elemInserted,edgeIndexInserted]=relabelElement(elem,edgeIndex,indexForSplit,node,edgeSet)

fprintf('Relabeling elements with labeled edges larger than 2...\n')
faceorder=[1 2 3; 1 2 4; 1 3 4; 2 3 4];

% get elements that need to be split
elemSplit=elem(indexForSplit,:);

% get centroids and add them into nodes
centroid=meshcentroid(node,elemSplit);
offset=size(node,1);
node=[node; centroid];

%% get new elements and new labelings

% 4 new elements will be generated for each element to be split
splitElem=cell(length(indexForSplit),1);  % store new elements after spliting
splitLocalEdge=cell(length(indexForSplit),1);  % store local edge labeling for the new elements
for i=1:size(elemSplit,1)
    % one element can be split into 4 new elements
    newElem=zeros(size(faceorder,1),4);  % 4 nodes in an element
    newLocalEdge=zeros(size(faceorder,1),6);  % 6 edges in an element
    for j=1:size(faceorder,1)
        % offset+i is the node index of the centroid
        newElem(j,:)=[offset+i,...
                        elemSplit(i,faceorder(j,1)),...
                        elemSplit(i,faceorder(j,2)),...
                        elemSplit(i,faceorder(j,3))];
        newLocalEdge(j,:)=labelSingleElem(newElem(j,:),edgeSet);
        if sum(logical(newLocalEdge(j,:)))>2
           warning('Split failed')
        end
    end
    splitElem{i}=newElem;
    splitLocalEdge{i}=newLocalEdge;
end

%% Insert new elements to original elements

if indexForSplit(1)==1
    elemInserted=[];
    edgeIndexInserted=[];
else
    elemInserted=elem(1:indexForSplit(1)-1,:);
    edgeIndexInserted=edgeIndex(1:indexForSplit(1)-1,:);
end
elemInserted=[elemInserted; splitElem{1}];
edgeIndexInserted=[edgeIndexInserted; splitLocalEdge{1}];

for i=2:length(indexForSplit)
    elemInserted=[elemInserted; elem(indexForSplit(i-1)+1:indexForSplit(i)-1,:)];
    elemInserted=[elemInserted; splitElem{i}];
    
    edgeIndexInserted=[edgeIndexInserted; edgeIndex(indexForSplit(i-1)+1:indexForSplit(i)-1,:)];
    edgeIndexInserted=[edgeIndexInserted; splitLocalEdge{i}];
end

if indexForSplit(end)<size(elem,1)
    elemInserted=[elemInserted; elem(indexForSplit(end)+1:end,:)];
    edgeIndexInserted=[edgeIndexInserted; edgeIndex(indexForSplit(end)+1:end,:)];
end
fprintf('Finished relabeling.\n')

%% function: label local edges for single element
function localEdge=labelSingleElem(elem,edgeSet)
% label a single element
n2e=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
localEdge=zeros(1,size(n2e,1));
c=1;
for i=1:size(n2e,1)
    % check if the edge is in the set
    curEdge=sort([elem(n2e(i,1)) elem(n2e(i,2))],2);
    if edgeSet.contains(num2str(curEdge))
        localEdge(c)=i;
        c=c+1;
    end
end


%% function: obtain the radii for the labeled edges
function radii=getVesselRadii(elem,index,noder)
n2e=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];  % node to edge mapping
radii=zeros(size(index,1),2);
for i=1:size(index,1)
    % edge 1
    if index(i,1)~=0
        ei1=n2e(index(i,1),:);
        radii(i,1)=max([noder(elem(i,ei1(1))),noder(elem(i,ei1(2)))]);
    end
    
    % edge 2
    if index(i,2)~=0
        ei2=n2e(index(i,2),:);
        radii(i,2)=max([noder(elem(i,ei2(1))),noder(elem(i,ei2(2)))]);
    end
end


%% function: build cfg for immc
function cfg=inputMesh2Cfg(elem,index,radii,node,noder)

% turn on iMMC
cfg.implicit=1;

% post-processing to fit to iMMC input format
index=index(:,1:4);
index(index==0)=7;
% 0 to 5 represet local edge labeled as vessel, 6 represent not labeled as
% vessel
index=index-1;

node=[node noder];
elem=[elem index radii zeros(size(radii,1),2)];

cfg.node=node;
cfg.elem=elem;
