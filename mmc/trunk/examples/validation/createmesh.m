sessionid='cube';

addpath('../../matlab/');
[node,elem]=genT5mesh(0:2:60,0:2:60,0:2:60);

vol=elemvolume(elem,node,'signed');
idx=find(vol<0);
elem(idx,[1 2 3 4])=elem(idx,[1 2 4 3]);

savemmcmesh(sessionid,node,elem,[]);
