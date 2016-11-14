sessionid='cube';

addpath('../../matlab/');
[node,face,elem]=meshabox([0 0 0],[60 60 20],2,2);
elem(:,1:4)=meshreorient(node,elem(:,1:4));

srcpos=[30.1,30.1,0];
savemmcmesh(sessionid,node,elem);
plotmesh(node,elem)
eid=tsearchn(node,elem(:,1:4),srcpos)
