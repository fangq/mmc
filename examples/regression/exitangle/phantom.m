clear;
phantom_l=50;
phantom_w=40;
phantom_h=21;
bar_dis=16.8;
bar_h=11;
bar_c1=[(phantom_l-bar_dis)/2,bar_h];
bar_c2=[(phantom_l+bar_dis)/2,bar_h];

%% old mesh

[node,face,elem] = meshabox([0 0 0],[phantom_l phantom_w phantom_h],1.2,10); 
elem(:,1:4) = meshreorient(node(:,1:3),elem(:,1:4));
elem(:,5) = 1;
plotmesh(node,elem);

% savemmcmesh('tank',node,elem);

%% add source

source = [5 5 23];
source_dir = [0 0 -1];
param1 = [40 0 0 0];
param2 = [0 30 0 0];
cfg=struct('srctype','planar','srcpos',source,'srcdir',source_dir,'srcparam1',param1,'srcparam2',param2);
source_domain = mmcsrcdomain(cfg,[min(node);max(node)]);

[newnode,newelem] = mmcaddsrc(node,elem,source_domain);
plotmesh(newnode,newelem);

savemmcmesh('tank',newnode,newelem);
