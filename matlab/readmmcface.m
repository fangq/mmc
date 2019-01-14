function face=readmmcface(filename)
%
% face=readmmcface(filename)
%
% Loading MMC surface triangle file
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     filename: the file name to the surface element data file
%
% output:
%     face: the surface triangle element list 
%
% example:
%     face=readmmcface('face_sph1.dat');
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

fid=fopen(filename,'rt');
[header,c]=fscanf(fid,'%d',2);
face=fscanf(fid,'%d',[5 header(2)]);
fclose(fid);

face=face';
face=face(:,2:4);
