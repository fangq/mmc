addpath ../../matlab
srcpos=[50.5 50.5 0];
fprintf('Generating mesh...');
[node,elem] = genT5mesh(0:4:100,0:4:100,0:4:80);
fprintf('done\nSaving mesh..');
savemmcmesh('dcs',node(:,1:3),elem); % 4th column of node is label and confuses savemmcmesh volume calculation
fprintf('done\nFinding element enclosing source: ');
eid=tsearchn(node(:,1:3),elem(:,1:4),srcpos);
fprintf('%d\n',eid);

