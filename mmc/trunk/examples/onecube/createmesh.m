sessionid='onecube';

addpath('../../matlab/');
[node,elem]=genT6mesh(0:1,0:1,0:1);
elem=sortrows(elem);
elem(:,1:4)=meshreorient(node,elem(:,1:4));

node=node*10;

srcpos=[2,8,0];
savemmcmesh(sessionid,node,elem,[]);
eid=tsearchn(node,elem(:,1:4),srcpos)

