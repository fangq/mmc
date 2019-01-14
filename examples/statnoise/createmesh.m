sessionid='cube20';

addpath('../../matlab/');
[node,elem]=genT5mesh(0:1:20,0:1:20,0:1:20);

elem(:,1:4)=meshreorient(node,elem(:,1:4));

srcpos=[10.1,10.2,0];
savemmcmesh(sessionid,node,elem,[]);
eid=tsearchn(node,elem(:,1:4),srcpos)

