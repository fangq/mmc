
% you may also use meshgrid5 in iso2mesh 1.7.9 or newer
% it is identical to genT5mesh

srcpos=[30.1 30.2 0.0];
[node,elem]=genT6mesh(0:2:60,0:2:60,0:2:60);
e0=tsearchn(node,elem,srcpos)
savemmcmesh('replaytest',node,elem);
