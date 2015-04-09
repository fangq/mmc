function [newnode,newelem]=mmcaddsrc(node,elem,src,varargin)
%
% [newnode,newelem]=mmcaddsrc(node,elem,src,opt)
%   or
% [newnode,newelem]=mmcaddsrc(node,elem,src,'Param1',value1, 'Param2',value2, ...)
%
% Adding an external wide-field (polyhedral) source domain to an
% existing tetrahedral mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%     node: the node list of the original mesh
%     elem: the elem list of the original mesh
%
% output:
%     newnode: the node list of the modified mesh, by default, 
%              newnode is a concatination of src at the end of node
%     newelem: the elem list of the modified mesh, by default, the elements
%              between the convex hull of the combined node/src and
%              the original mesh surface are labeled as 0; the elements
%              sharing the source domain are labeled as -1.
%     opt: optional, a struct defining additional options, the possible fields
%              are identical to those supported by meshrefine.m in iso2mesh.
%              the opt struct can be replaced by a list of 'param','value' pairs
% example:
%     [node,face,elem]=meshasphere([0 0 0],24,1,2);
%     elem(:,5)=1;
%
%     source=double([-5 -5;-5 5;5 5;5 -5]);
%     source(:,3)=30;
%     [newnode,newelem]=mmcaddsrc(node,elem,source);
%     plotmesh(newnode,newelem,'x>-5');
%
%     % example with additional options
%     [newnode,newelem]=mmcaddsrc(node,elem,source,'extcmdopt','-Y -a1','extlabel',-2);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

opt=varargin2struct(varargin{:});
if(~isfield(opt,'extcmdopt'))
   opt.extcmdopt='-Y';
end
opt.newnode=src;
[newnode,newelem]=meshrefine(node,elem,opt);
