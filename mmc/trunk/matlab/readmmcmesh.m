function varargout=readmmcmesh(key)
%
% [node elem]=readmmcmesh(key)
%
% Loading MMC node and element data files
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     key: the file name stub to the node coordinate file. The full 
%          file names are {node,elem}_key.dat
%
% output:
%     node: the node coordinate list
%     elem: the tetrahedra node index list
%
% example:
%     [node elem]=readmmcmesh('sph1');
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if(nargout>=1)
    varargout{1}=readmmcnode(['node_',key,'.dat']);
end
if(nargout>=2)
    varargout{2}=readmmcelem(['elem_',key,'.dat']);
end
if(nargout>=3)
    varargout{3}=readmmcelem(['elem_',key,'.dat']);
end
