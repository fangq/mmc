sessionid='onecube';

addpath('../../matlab/');
[node,elem]=genT6mesh(0:1,0:1,0:1);
elem=sortrows(elem);
elem(:,1:4)=meshreorient(node,elem(:,1:4));
elem(:,5)=1;

node=node*10;

srcpos=[2,8,0];
eid=tsearchn(node,elem(:,1:4),srcpos)
elem(eid,5)=2;
savemmcmesh(sessionid,node,elem,[]);

