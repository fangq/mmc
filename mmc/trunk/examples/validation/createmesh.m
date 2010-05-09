sessionid='cube';

addpath('../../matlab/');
[node,elem]=genT5mesh(0:1:60,0:1:60,0:1:60);

elem(:,1:4)=meshreorient(node,elem(:,1:4));

srcpos=[30.1,30.2,0];
savemmcmesh(sessionid,node,elem,[]);
eid=tsearchn(node,elem(:,1:4),srcpos)

