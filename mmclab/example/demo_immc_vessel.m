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
    % vnode, velem: nodes and elements for generated mesh
    [vnode, velem]=s2m(node,num2cell(face,2),1,100,'tetgen1.5',[],[],'-YY');

    %% 2. Create edgeroi
    % (a) get mapping from edge to radius (you can build your own e2r)
    e2r=edge2radius(edge,noder);
    % (b) get edgeroi
    [vnode,velem,eroi]=edgeroi(vnode,velem,e2r);

    %% 3. Run mmc

    res=0.001;  % resolution: 0.001mm
    cfg.nphoton=1e6;
    cfg.srctype='pencil';
    cfg.srcpos=[150.5 150.5 605]*res;
    cfg.srcdir=[0 0 -1];
    cfg.elemprop=ones(size(velem,1),1);
    cfg.prop=[0 0 1 1;
        0.0458   35.6541    0.9000    1.3700;   % dermis
        23.0543    9.3985    0.9000    1.3700]; % blood
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-9;
    cfg.debuglevel='TP';
    cfg.method='grid';
    cfg.steps=[0.001 0.001 0.001];
    cfg.gpuid=-1;
    cfg.isreflect=1;
    cfg.outputtype='energy';
    
    cfg.elem=velem;
    cfg.node=vnode*res;
    cfg.edgeroi=eroi*res;
    cfg.noderoi=zeros(size(vnode,1),1);
    cfg.noderoi(1:size(noder,1))=noder*res;
end

tic
flux=mmclab(cfg);
toc

%% plotting
figure,imagesc(log10(rot90(squeeze(flux.data(:,150,:)))))
axis equal;
title('energy deposition map');


%% function: get vessel radii
function e2r = edge2radius(edge,noder)

% vr = vesselradii(edge,noder)
% 
% get the radii for the input edges (vessels)
% 
% input:
%     edge: ne x 2 array. ne is the total number of edges. 2 columns represent
%     the indices of nodes connecting the edges
%     noder: nr x 1 array storing the radius of each node. nr is the total 
%     number of nodes
% output:
%     vr: ne x 1 array storing the radius of each edge

e2r=java.util.HashMap;  % set storing all vessel edges
edge=sort(edge,2);
for i=1:size(edge,1)
    vr=max([noder(edge(i,1)),noder(edge(i,2))]);
    e2r.put(num2str(edge(i,:)),vr);
end

%% function: get edgeroi
function [node,elem,eroi]=edgeroi(node,elem,e2r)
% 
% build the edgeroi for iMMC
% 
% input:
%     node: node coordinates for elem
%     elem: input mesh containing the edges
%     e2r: java.util.HashMap mapping edge to its radius
%          <key, value> -> <edge, radius>
%          Example:
%          edge=[1 2], radius=0.1
%          To build the mapping, you should do
%          e2r.put(num2str(edge),radius)
% output:
%     eroi: ne x 6 array, storing the cylindrial edge-roi radii in the 
%           nchoosek order [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]. 
%           A 0 indicates no ROI. Up to 2 edges are supported

eo=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % edge order in an element
eroi = zeros(size(elem,1),6);

i=1;
while i<=size(elem,1)
    enum = 0; % number of edges being labeled
    for j=1:6
        edge=sort([elem(i,eo(j,1)) elem(i,eo(j,2))]);
        if e2r.containsKey(num2str(edge))
            eroi(i,j)=e2r.get(num2str(edge));
            enum=enum+1;
        end
    end
    
    if enum>2
        % if there are more than 2 labeled edges,split and relabel the 
        % current element
        [node,elem,eroi]=splitelem(node,elem,eroi,e2r,i);
        i=i+3; % replace the old element with four new elements
    end
    
    i=i+1;
end

%% function
function [node,elem,eroi]=splitelem(node,elem,eroi,e2r,index)
% 
% split the element that has over 2 labeled edges
% 
% input:
%     node:  node coordinates for elem
%     elem:  input mesh containing the edges
%     eroi:  input data for cfg.edgeroi
%     e2r:   java.util.HashMap mapping edge to its radius
%     index: index of element that needs to be split
% 
% output:
%     elem:  mesh after spliting

eo=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % edge order in an element
fo=[1 2 3; 1 2 4; 1 3 4; 2 3 4]; % face order in an element

selem=elem(index,:); % element to be split
c=meshcentroid(node,selem); % centroid of the element to be split
cindex=size(node,1)+1;
node=[node; c]; % put the centroid in the end

newelem=[]; % 4 elements that have been split using centroid
newroi=[];  % new edge roi
for i=1:size(fo,1)
    e=[selem(fo(i,1)) selem(fo(i,2)) selem(fo(i,3)) cindex];
    newelem=[newelem; e];
    
    r=zeros(1,size(eo,1));
    % check the edge that has radius and fill it in edgeroi
    for j=1:size(eo,1)
        edge=sort([e(eo(j,1)),e(eo(j,2))]);
        if e2r.containsKey(num2str(edge))
            r(j)=e2r.get(num2str(edge));
        end
    end
    newroi=[newroi; r];
end

% remove the old element and insert the new elements
elem=[elem(1:index-1,:); newelem; elem(index+1:end,:)];
% remove the old roi and insert the new rois
eroi=[eroi(1:index-1,:); newroi; eroi(index+1:end,:)];