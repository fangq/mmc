addpath('../../matlab/');
node=readmmcnode('node_onecube.dat');
elem=readmmcelem('elem_onecube.dat');

load ad.txt
load mov.txt

figure;
hh=tetramesh(elem,node);
set(hh,'facealpha',0.1)
hold on
photonnum=max(mov(:,end-1));
for i=0:photonnum
     idx=find(mov(:,end-1)==i);
     plot3(mov(idx,2),mov(idx,3),mov(idx,4),'.-');
end
idx=find(mov(:,1)==0);
plot3(mov(idx,2),mov(idx,3),mov(idx,4),'co');

idx=find(mov(:,1)==2);
plot3(mov(idx,2),mov(idx,3),mov(idx,4),'+');

idx=find(mov(:,1)==3);
plot3(mov(idx,2),mov(idx,3),mov(idx,4),'co');

plot3(ad(:,1),ad(:,2),ad(:,3),'r.');
