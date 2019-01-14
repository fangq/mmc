no=readmmcnode('node_sphere.dat');
el=readmmcelem('elem_sphere.dat');
load planar.dat;
vv=reshape(planar(:,2),size(no,1),size(planar,1)/size(no,1));
vvv=sum(vv,2);
plotmesh([no,vvv],el(el(:,end)>0,1:4))
ppl=loadmch('planar.mch');
hold on
plotmesh(ppl(:,4:6),'r.'); % plot detected photon locations from a widefield detector
figure
plotmesh([no,log10(vvv)],el(el(:,end)>0,1:4),'x=0');
shading interp

