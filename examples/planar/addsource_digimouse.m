if(~exist('Digimouse_Mesh_1L.mat','file'))
     if(~exist('digimouse_mesh_version_1L.tar.gz','file'))
         websave('digimouse_mesh_version_1L.tar.gz','http://downloads.sourceforge.net/project/mcx/mmc/Digimouse%20FEM%20Mesh/Version%201/digimouse_mesh_version_1L.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fmcx%2Ffiles%2Fmmc%2FDigimouse%2520FEM%2520Mesh%2FVersion%25201%2F&ts=1450931781&use_mirror=skylineservers');
     end
     if(~exist('digimouse_mesh_version_1L.tar','file'))
         gunzip('digimouse_mesh_version_1L.tar.gz');
     end
     if(~exist('digimouse_mesh/Digimouse_Mesh_1L.mat','file'))
         untar('digimouse_mesh_version_1L.tar');
         movefile('digimouse_mesh/Digimouse_Mesh_1L.mat','./Digimouse_Mesh_1L.mat');
     end
end
load Digimouse_Mesh_1L;
elem(:,5)=1;

cfg=struct('srctype','planar','srcpos',[15 50 25],'srcdir',[0 0 -1],...
           'srcparam1',[10 0 0 0],'srcparam2',[0 10 0 0]);

[newnode,newelem]=mmcaddsrc(node,elem,cfg);
savemmcmesh('digimouse',newnode,newelem);
