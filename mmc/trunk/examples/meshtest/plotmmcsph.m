node1=readmmcnode('node_sph1.dat');  
elem1=readmmcelem('elem_sph1.dat');
load sph1.dat
sph1=reshape(sph1(:,end),[size(node1,1),50]);
s1=sum(sph1,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem1(:,1:4),node1,s1,[0 30.2 0; 0 30.2 1; 1 30.2 0]);	  
patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
view([0 1 0])
axis equal   

node2=readmmcnode('node_sph2.dat');						    
elem2=readmmcelem('elem_sph2.dat');										    
load sph2.dat
sph2=reshape(sph2(:,end),[size(node2,1),50]);
s2=sum(sph2,2);
figure
[cutpos,cutvalue,facedata]=qmeshcut(elem2(:,1:4),node2,s2,[0 30.2 0; 0 30.2 1; 1 30.2 0]);
patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
view([0 1 0])

node3=readmmcnode('node_sph3.dat');
elem3=readmmcelem('elem_sph3.dat');
load sph3.dat
sph3=reshape(sph3(:,end),[size(node3,1),50]);
s3=sum(sph3,2);
figure
[cutpos,cutvalue,facedata]=qmeshcut(elem3(:,1:4),node3,s3,[0 30.2 0; 0 30.2 1; 1 30.2 0]);
patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
view([0 1 0])
