function elem=readmmcelem(filename)
%
% elem=readmmcelem(filename)
%
% Loading MMC mesh element file
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     filename: the file name to the element data file
%
% output:
%     elem: the tetrahedral mesh element list 
%
% example:
%     elem=readmmcelem('elem_sph1.dat');
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

fid=fopen(filename,'rt');
[header,c]=fscanf(fid,'%d',2);
elem=fscanf(fid,'%d',[6 header(2)]);
fclose(fid);

elem=elem';
elem=elem(:,2:end);
