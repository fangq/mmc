figure
load enterpos.txt
subplot(121);
plotmesh(enterpos(:,1:3),'.')
load exitpos.txt 
subplot(122);
plotmesh(exitpos(:,1:3),'.') 

