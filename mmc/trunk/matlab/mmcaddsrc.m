function [newnode,newelem]=mmcaddsrc(node,elem,src,varargin)
%
% [newnode,newelem]=mmcaddsrc(node,elem,src,opt)
%   or
% [newnode,newelem]=mmcaddsrc(node,elem,src,'Param1',value1, 'Param2',value2, ...)
%
% Adding an external wide-field (polyhedral) source domain to an
% existing tetrahedral mesh
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     node: the node list of the original mesh
%     elem: the elem list of the original mesh
%     src:  the coordinates of the vertices that define the source domain
%           if src is a struct, it will be treated as an MMCLAB input
%           structure and the source vertices will be computed
%           automatically
%     opt: optional, a struct defining additional options, the possible fields
%              are identical to those supported by meshrefine.m in iso2mesh.
%              the opt struct can be replaced by a list of 'param','value' pairs
%
% output:
%     newnode: the node list of the modified mesh, by default, 
%              newnode is a concatination of src at the end of node
%     newelem: the elem list of the modified mesh, by default, the elements
%              between the convex hull of the combined node/src and
%              the original mesh surface are labeled as 0; the elements
%              sharing the source domain are labeled as -1.
%
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
%     % example with mmclab cfg input
%     cfg=struct('srctype','cone','srcpos',[0 0 28],'srcdir',[0 0 -1]);
%     [newnode,newelem]=mmcaddsrc(node,elem,cfg);
%     figure; plotmesh(newnode,newelem);
%
%     % for box-shaped objects, you have to use meshabox with matlab 7.11
%     % or newer, older matlab versions will fail to get the convex hull
%     [node,face,elem]=meshabox([-10 -10 -10],[10 10 0],1,2);
%     elem(:,5)=1;
%
%     cfg=struct('srctype','fourier','srcpos',[-5 -5 10],'srcdir',...
%        [0 0 -1],'srcparam1',[10 0 0 3],'srcparam2',[0 8 0 1]);
%     [newnode,newelem]=mmcaddsrc(node,elem,cfg);
%     figure; plotmesh(newnode,newelem);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

opt=varargin2struct(varargin{:});
if(~isfield(opt,'extcmdopt'))
   opt.extcmdopt='-Y';
end
if(isstruct(src))
    src=mmcsrcdomain(src,[min(node);max(node)],opt);
end
opt.newnode=src;
[newnode,newelem]=meshrefine(node,elem,opt);
