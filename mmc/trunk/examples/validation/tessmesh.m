sessionid='cube';

addpath('../../matlab/');
%[node,elem]=genT5mesh(0:2:60,0:2:60,0:2:60);

[no,fc]=meshabox([0 0 0], [60 60 60], 100);
[no,fc]=removeisolatednode(no,fc);

ISO2MESH_TETGENOPT='-YY'

[node,elem]=s2m(no,fc,1,100);

elem(:,1:4)=meshreorient(node,elem(:,1:4));

srcpos=[30.1,30.2,0];
savemmcmesh(sessionid,node,elem,[]);
eid=tsearchn(node,elem(:,1:4),srcpos)

