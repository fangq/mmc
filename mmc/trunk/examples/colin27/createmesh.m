addpath('../../../matlab/'); 

if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.tar.gz','file'))
          urlwrite('http://downloads.sourceforge.net/project/mcx/mmc/AdultBrain%20Atlas%20FEM%20Mesh/Version%202/MMC_Collins_Atlas_Mesh_Version_2L.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fmcx%2Ffiles%2Fmmc%2FAdultBrain%2520Atlas%2520FEM%2520Mesh%2FVersion%25202%2F&ts=1451344236&use_mirror=tcpdiag','MMC_Collins_Atlas_Mesh_Version_2L.tar.gz');
     end
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.tar','file'))
          gunzip('MMC_Collins_Atlas_Mesh_Version_2L.tar.gz');
     end
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
          untar('MMC_Collins_Atlas_Mesh_Version_2L.tar');
     end
end

load MMC_Collins_Atlas_Mesh_Version_2L.mat
savemmcmesh('brain',node,elem);
