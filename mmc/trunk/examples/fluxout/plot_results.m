
mua=0.01;
tstep=1e-10;
gatenum=50;

node=readmmcnode('node_cube.dat');
elem=readmmcelem('elem_cube.dat');
face=readmmcface('face_cube.dat');

data=load('cube.dat');
data=reshape(data(:,2),[],gatenum);
cwdata=sum(data,2);
Eabsorb=sum(cwdata.*nodevolume(node,elem))*tstep*mua;

facedata=load('cube_face.dat');
facedata=facedata(:,2);
facedata=reshape(facedata,[],gatenum);
cwfacedata=sum(facedata,2);
Eout=sum(cwfacedata)*tstep;

Etotal=Eabsorb+Eout;

facelabel=round(log10(cwfacedata+eps));

figure(1);
plotmesh(node,[face,facelabel])

figure(2);
plotmesh(node,log10(cwfacedata+eps));
