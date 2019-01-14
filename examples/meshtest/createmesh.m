%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Create meshes for a sphere inside a cubic domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation

% 1. you need to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

% 2. you need to add the path to MMC matlab folder
addpath('../../matlab')


% create a surface mesh for a 10mm radius sphere
[no,el]=meshasphere([30 30 30],10,1.0);


% generate a coarse volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 20

ISO2MESH_SESSION='mmcmesh2_';

srcpos=[30. 30. 0.];
fixednodes=[30.,30.,0.05; 30 30 30];
nodesize=[ones(size(no,1),1) ; 0.5; 3];
nfull=[no;fixednodes];
[node3,elem3,face3]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,8,[30 30 30],[],[1.5 1.5 1.5 1.5 5 5 5 5]);
[node3,elem3]=sortmesh(srcpos,node3,elem3,1:4);
elem3(:,1:4)=meshreorient(node3,elem3(:,1:4));
elem3(:,5)=elem3(:,5)+1;
savemmcmesh('mesh2',node3,elem3(:,1:5),[]);
eid3=tsearchn(node3,elem3(:,1:4),srcpos);

% generate a dense volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 5

ISO2MESH_SESSION='mmcmesh1_';

nodesize=[1*ones(size(no,1),1) ; 1; 1];
[node2,elem2,face2]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,2,[30 30 30],[],[1 1 1 1 1 1 1 1]);
[node2,elem2]=sortmesh(srcpos,node2,elem2,1:4);
elem2(:,1:4)=meshreorient(node2,elem2(:,1:4));
elem2(:,5)=elem2(:,5)+1;
savemmcmesh('mesh1',node2,elem2(:,1:5),[]);
eid2=tsearchn(node2,elem2(:,1:4),srcpos);

% reduce the surface node numbers to 20%

ISO2MESH_SESSION='mmcmesh0_';

%[no2,el2]=meshresample(no,el,0.2);

% using the coarse spherical surface, we generate a coarse volumetric
% mesh with maximum volume of 10

nodesize=[ones(size(no,1),1); 3];
nfull=[no; 30 30 30];
[node1,elem1,face1]=surf2mesh([nfull,nodesize],el,[0 0 0],[60.1 60.1 60.1],1,10,[30 30 30],[],[2 2 2 2 5 5 5 5]);
[node1,elem1]=sortmesh(srcpos,node1,elem1,1:4);
elem1(:,1:4)=meshreorient(node1,elem1(:,1:4));
elem1(:,5)=elem1(:,5)+1;
savemmcmesh('mesh0',node1,elem1(:,1:5),[]);
eid1=tsearchn(node1,elem1(:,1:4),srcpos);

clear ISO2MESH_SESSION

fid=fopen('initial_elem.txt','wt');
fprintf(fid,'mesh0: %d\nmesh1: %d\nmesh2: %d\n',eid1,eid2,eid3);
fclose(fid);
