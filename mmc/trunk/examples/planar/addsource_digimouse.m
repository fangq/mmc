load Digimouse_Mesh_1L;
elem(:,5)=1;

cfg=struct('srctype','planar','srcpos',[15 50 25],'srcdir',[0 0 -1],...
           'srcparam1',[10 0 0 0],'srcparam2',[0 10 0 0]);

[newnode,newelem]=mmcaddsrc(node,elem,cfg);
savemmcmesh('digimouse',newnode,newelem);
