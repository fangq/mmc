function varargout=mmcadddet(varargin)
%
% [newnode,newelem]=mmcadddet(node,elem,det,opt)
%   or
% [newnode,newelem]=mmcadddet(node,elem,det,'Param1',value1, 'Param2',value2, ...)
%
% Adding an external wide-field (polyhedral) detector domain to an
% existing tetrahedral mesh
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     see mmcaddsrc for details.
%
% output:
%     see mmcaddsrc for details.
%
% example:
%     [node,face,elem]=meshasphere([0 0 0],24,1,2);
%     elem(:,5)=1;
%
%     cfg=struct('srctype','cone','srcpos',[0 0 28],'srcdir',[0 0 -1]);
%     [nodesrc,elemsrc]=mmcaddsrc(node,elem,cfg);
%
%     % example with additional options
%     cfg=struct('srctype','planar','srcpos',[-30 -5 5],'srcdir',[1 0 0],...
%                'srcparam1',[0 10 0 0],'srcparam2',[0 0 8 0]);
%     [nodedet,elemdet]=mmcadddet(nodesrc,elemsrc,cfg);
%     plotmesh(nodedet,elemdet);
%
%     cfg=struct('srctype','disk','srcpos',[-30 0 0],'srcdir',[1 0 0],'srcparam1',[4 0 0 0]);
%     [nodedet,elemdet]=mmcadddet(nodesrc,elemsrc,cfg);
%     figure;plotmesh(nodedet,elemdet);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

[varargout{1:nargout}]=mmcaddsrc(varargin{:},'extcorelabel',-2,'KeepShape',1,'Expansion',1.0);

