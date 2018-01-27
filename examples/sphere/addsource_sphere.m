[node,face,elem]=meshasphere([0 0 0],24,1,2);
elem(:,5)=1;

cfg=struct('srctype','planar','srcpos',[-5 -5 25],'srcdir',[0 0 -1],...
           'srcparam1',[10 0 0 0],'srcparam2',[0 10 0 0]);
[newnode,newelem]=mmcaddsrc(node,elem,cfg);

cfg=struct('srctype','planar','srcpos',[-30 -5 -5],'srcdir',[1 0 0],...
           'srcparam1',[0 10 0 0],'srcparam2',[0 0 10 0]);
[nodedet,elemdet]=mmcadddet(newnode,newelem,cfg);
plotmesh(nodedet,elemdet);

savemmcmesh('sphere',nodedet,elemdet);
