
%% original mesh

[node,face,elem]=meshabox([0 0 0],[60 60 20],2, 2);
elem(:,5)=1;

figure(1);
plotmesh(node,elem);

%% extended source domain

srcdef=struct('srctype','planar','srcpos',[10,10,-2],'srcdir',[0 0 1],...
                'srcparam1',[40 0 0 40],'srcparam2',[0 40 0 40]);
            
[newnode,newelem]=mmcaddsrc(node,elem,srcdef);

figure(2);
plotmesh(newnode,newelem);

%% extended detector domain

detdef=struct('srctype','planar','srcpos',[10,10,20.1],'srcdir',[0 0 -1],...
                'srcparam1',[40 0 0 40],'srcparam2',[0 40 0 40]);
            
[newnode2,newelem2]=mmcadddet(newnode,newelem,detdef);

figure(3);
plotmesh(newnode2,newelem2);

savemmcmesh('replaywide',newnode2,newelem2);