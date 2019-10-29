addpath('../../../matlab/'); 

if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
     f=websave('MMC_Collins_Atlas_Mesh_Version_2L.mat','https://github.com/fangq/Colin27BrainMesh/raw/master/MMC_Collins_Atlas_Mesh_Version_2L.mat');
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
          warning('Please download Colin27 atlas mesh from https://github.com/fangq/Colin27BrainMesh/raw/master/MMC_Collins_Atlas_Mesh_Version_2L.mat')
     end
end

load MMC_Collins_Atlas_Mesh_Version_2L.mat
savemmcmesh('brain',node,elem);
