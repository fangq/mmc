%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Create meshes for a sphere inside a cubic domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation

% 1. you have to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

% 2. you have to add the path to MMC matlab folder
addpath('../../matlab')


% create a gray-scale field
dim=60;
[xi,yi,zi]=meshgrid(0:0.5:dim,0:0.5:dim,0:0.5:dim);
dist=sqrt((xi-30).^2+(yi-30).^2+(zi-30).^2);
clear xi yi zi;

% extract a level-set at v=20, being a sphere with R=20
% the maximum element size of the surface triangles is 2

[v0,f0]=vol2restrictedtri(dist,10,[60 60 60],100*20,30,2,2,40000);
v0=(v0-0.5)*0.5;

% iso2mesh will also produce a surface for the bounding box, remove it
facecell=finddisconnsurf(f0);
sphsurf=facecell{1};

if( sum((v0(sphsurf(1,1),:)-[30 30 30]).^2) > 25*25)
   sphsurf=facecell{2};
end
%trisurf(sphsurf,v0(:,1),v0(:,2),v0(:,3));
%axis equal;
idx=unique(sphsurf);  % this is the index of all the nodes on the sphere

% show the histogram of the displacement error for the nodes on the sphere
r0=sqrt((v0(idx,1)-30).^2+(v0(idx,2)-30).^2+(v0(idx,3)-30).^2);
%figure;hist(r0,100);

% we only take the nodes on the surface
[no,el]=removeisolatednode(v0,sphsurf);
[no,el]=meshcheckrepair(no,el);

% generate a coarse volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 20

ISO2MESH_SESSION='mmcsph3_';

srcpos=[30. 30. 0.];
fixednodes=[30.,30.,0.05; 30 30 30];
nodesize=[ones(size(no,1),1) ; 0.3; 3];
nfull=[no;fixednodes];
[node3,elem3,face3]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,8,[30 30 30],[],[1.5 1.5 1.5 1.5 6 6 6 6]);
[node3,elem3]=sortmesh(srcpos,node3,elem3,1:4);
elem3(:,1:4)=meshreorient(node3,elem3(:,1:4));
savemmcmesh('sph3',node3,elem3(:,1:5),[]);
eid3=tsearchn(node3,elem3(:,1:4),srcpos);

% generate a dense volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 5

ISO2MESH_SESSION='mmcsph2_';

nodesize=[0.7*ones(size(no,1),1) ; 0.2; 2];
[node2,elem2,face2]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,2,[30 30 30],[],[1 1 1 1 5 5 5 5]);
[node2,elem2]=sortmesh(srcpos,node2,elem2,1:4);
elem2(:,1:4)=meshreorient(node2,elem2(:,1:4));
savemmcmesh('sph2',node2,elem2(:,1:5),[]);
eid2=tsearchn(node2,elem2(:,1:4),srcpos);

% reduce the surface node numbers to 20%

ISO2MESH_SESSION='mmcsph1_';

[no2,el2]=meshresample(no,el,0.2);

% using the coarse spherical surface, we generate a coarse volumetric
% mesh with maximum volume of 10

nodesize=[2*ones(size(no2,1),1)];
[node1,elem1,face1]=surf2mesh([no2,nodesize],el2,[0 0 0],[60.1 60.1 60.1],1,10,[30 30 30],[],[1 1 1 1 5 5 5 5]);
[node1,elem1]=sortmesh(srcpos,node1,elem1,1:4);
elem1(:,1:4)=meshreorient(node1,elem1(:,1:4));
savemmcmesh('sph1',node1,elem1(:,1:5),[]);
eid1=tsearchn(node1,elem1(:,1:4),srcpos);

clear ISO2MESH_SESSION

fid=fopen('initial_elem.txt','wt');
fprintf(fid,'sph1: %d\nsph2: %d\nsph3: %d\n',eid1,eid2,eid3);
fclose(fid);
