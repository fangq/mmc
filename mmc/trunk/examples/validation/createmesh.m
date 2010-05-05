sessionid='cube';

addpath('../../matlab/');
[node,elem]=genT5mesh(0:2:60,0:2:60,0:2:60);

elem(idx,1:4)=reorient(node,elem(:,1:4));

savemmcmesh(sessionid,node,elem,[]);
