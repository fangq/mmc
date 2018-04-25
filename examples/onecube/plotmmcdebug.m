addpath('../../matlab/');
node=readmmcnode('node_onecube.dat');
elem=readmmcelem('elem_onecube.dat');

load ad.txt
load mov.txt
srcpos=[2 8 0];

hh=tetramesh(elem(:,1:4),node);
set(hh,'facealpha',0.1)
hold on

photonnum=max(mov(:,end-1));
tracks=[];
for i=0:photonnum
     idx=find(mov(:,end-1)==i);
     tracks=[tracks; [srcpos(1);mov(idx,2)],[srcpos(2);mov(idx,3)],[srcpos(3);mov(idx,4)]; nan nan nan];
end
plot3(tracks(:,1),tracks(:,2),tracks(:,3),'.-');
idx=find(mov(:,1)==0);
plot3(mov(idx,2),mov(idx,3),mov(idx,4),'ko');

idx=find(mov(:,1)==2);
plot3(mov(idx,2),mov(idx,3),mov(idx,4),'+');

plot3(ad(:,1),ad(:,2),ad(:,3),'r.');
