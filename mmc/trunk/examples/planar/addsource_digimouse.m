clear;

load Digimouse_Mesh_1L;
elem(:,5)=1;

source=double([15 50 25;15 60 25;25 50 25;25 60 25]);

% new nodes are exactly four corners of the plane
newnode=source;

allnode=[node;newnode];

% save source node id of new nodes for future use
snid=[length(node)+1:length(allnode)];

% create the convex hull
outface=convhulln(allnode);
outface=sort(outface,2);

% the inner face of the domain to be meshed
face=volface(elem(:,1:4));
inface=sort(face(:,1:3),2);

% remove coinside triangles
bothsides=removedupelem([outface;inface]);

% define a seed point to avoid meshing the interior space
holelist=surfseeds(node,face(:,1:3));

% let tetgen not to insert point on the surface to be consistent
ISO2MESH_TETGENOPT='-YY';

% generate the tetrahedral mesh for the space between inface and outface
[no,el,fc]=surf2mesh(allnode,bothsides,[],[],1,10,[],holelist);
clear ISO2MESH_TETGENOPT

% create map between all nodes and newly generated nodes
[isinside,map]=ismember(no,allnode,'rows');

% if all the nodes are mapped successfully, change the node number of list el
if numel(find(isinside==0))==0
    el2=zeros(size(el));
    for i=1:length(no)
        el2(el==i)=map(i);
    end
    % label all new elements with -1
    el2(:,5)=0;
    
    % merge nodes/elements and replace the original ones
    node=allnode;
    elem=[elem;el2];
    
    % search elements that contain source(s) and save their id(s)
    comb=sort(nchoosek(snid,3),2);  % all possible combinations of source triangles
    elm=sort(elem(:,1:4),2);        % all sorted elements
    [iselm,seid]=ismember(comb,elm(:,2:end),'rows');
    seid=seid(iselm);               % source element id list
    elem(seid,5)=-1;
    savemmcmesh('digimouse',node,elem);
    display('All nodes mapped successfully');
else
    display('Fail to map all new-generated nodes');
end