
%% original mesh

[node,face,c0]=latticegrid([0 60],[0 60],[0 5 10]);
c0(:,4)=[2;3];   % maximum element size for bottom (label 1) and top (label 2) layers
[node,elem]=surf2mesh(node,face,[],[],1,[],c0);

figure(1);
plotmesh(node,elem);

srcpos = [25.0 35.0 10.0];

%% extended detector domain

detdef=struct('srctype','planar','srcpos',[10,10,-1],'srcdir',[0 0 1],...
                'srcparam1',[40 0 0 40],'srcparam2',[0 40 0 40]);
            
[newnode,newelem]=mmcadddet(node,elem,detdef);

figure(2);
plotmesh(newnode,newelem);

newelem(:,1:4)=meshreorient(newnode,newelem(:,1:4));
eid=tsearchn(newnode,newelem(:,1:4),srcpos)
savemmcmesh('widedet',newnode,newelem);