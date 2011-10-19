function [mua,mus,g,n]=load_mc_prop(fname)
%
%   [mua,mus,g,n]=load_mc_prop(fname)
%
%   Loads the absorption and scattering coefficients, as well as the 
%   scattering anisotropy and index of refraction from the MMC optical
%   property file given as an argument
%  
%   author: Stefan Carp (carp <at> nmr.mgh.harvard.edu


fid=fopen(fname,'r');
if fid<0, 
  error('File %s not found\n',fname);
end

mc_version=fscanf(fid,'%f',1); if mc_version>1, fprintf('WARNING: MMC version greater than 1\n'); end
num_tissues=fscanf(fid,'%f',1);

for I=1:num_tissues,
    tiss_ind=fscanf(fid,'%f',1);
    mua(tiss_ind)=fscanf(fid,'%f',1);
    mus(tiss_ind)=fscanf(fid,'%f',1);
    g(tiss_ind)=fscanf(fid,'%f',1);
    n(tiss_ind)=fscanf(fid,'%f',1);
end

