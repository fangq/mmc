function [srcnode,srcface]=mmcsrcdomain(cfg,meshbbx,varargin)
%
% [srcnode,srcface]=mmcsrcdomain(cfg)
%
% Defining a source domain (for launching new photons) in the form of 
% polyhedra based on an MMCLAB simulation configuration structure
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%     cfg: the simulation configuration structure used by mmclab, the
%          following fields are required
%               cfg.srcpos, cfg.srcdir, cfg.srctype, cfg.srcparam1,
%               cfg.srcparam2
%          run "help mmclab" for details.
%     meshbbx: the bounding box of the original mesh in the form of
%          [min(node); max(node)], i.e. a 2x3 matrix
%
% output:
%     srcnode: the vertices of the source domain
%     srcface (optional): the polyhedral patches of the source domain
%
% example:
%     [node,face,elem]=meshasphere([0 0 0],20,1);
%     elem(:,5)=1;
%     cfg=struct('srctype','cone','srcpos',[0 0 25],'srcdir',[0 0 -1]);
%     srccone=mmcsrcdomain(cfg,[min(node);max(node)]);
%
%     cfg=struct('srctype','fourier','srcpos',[-5 -5 25],'srcdir',...
%        [0 0 -1],'srcparam1',[10 0 0 3],'srcparam2',[0 8 0 1]);
%     srcfourier=mmcsrcdomain(cfg,[min(node);max(node)]);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

% the launching domain should be slightly larger than the actual source
% domain to avoid starting from an edge or vertex
expansion=1.1; 

opt=varargin2struct(varargin{:});
if(~isstruct(cfg))
   error('input cfg must be a struct');
end
if(~isfield(cfg,'srcpos') || ~isfield(cfg,'srcdir') || ...
   ~isfield(cfg,'srctype') )
    error('you must at least define cfg.srcpos, cfg.srcdir, cfg.srctype');
end

if(strcmp(cfg.srctype,'isotropic'))
    cfg.srcdir=rand(1,3);
    cfg.srcdir=cfg.srcdir/norm(cfg.srcdir);
end
cfg.srcdir=cfg.srcdir/norm(cfg.srcdir);
domainradius=min(meshbbx(2,:)-meshbbx(1,:))*0.25;
if(strcmp(cfg.srctype,'pencil'))
   warning(['for external pencil beams, you are highly recommended to set ' ...
       'cfg.e0=''>'', instead of using this script']); 
end

if(strcmp(cfg.srctype,'pencil') || strcmp(cfg.srctype,'cone') || ...
        strcmp(cfg.srctype,'zgaussian')|| strcmp(cfg.srctype,'isotropic') || strcmp(cfg.srctype,'arcsine'))
    srcnode=orthdisk(cfg.srcpos,cfg.srcpos+cfg.srcdir,domainradius,3);
    srcface=[1 2 3];
elseif(strcmp(cfg.srctype,'planar') || strcmp(cfg.srctype,'pattern') || ...
        strcmp(cfg.srctype,'fourier') ||strcmp(cfg.srctype,'fourierx')||strcmp(cfg.srctype,'fourier2d'))
    v0=cfg.srcpos(1:3);
    v1=cfg.srcparam1(1:3);
    v2=cfg.srcparam2(1:3);
    if(strcmp(cfg.srctype,'fourierx')||strcmp(cfg.srctype,'fourier2d'))
        v2=cross(v0,v1);
    end
    voff1=(v1+v2)*0.5;
    voff2=(v1-v2)*0.5;
    if(jsonopt('SplitDomain',0,opt))
        srcnode=[v0; v0+v1; v0+v1+v2; v0+v2]+[-voff1; voff2; voff1; -voff2]*(expansion-1.0);
        srcface=[1 2 3;3 4 1];
    else
        srcnode=orthdisk(cfg.srcpos+voff1,cfg.srcpos+voff1+cfg.srcdir,(max(norm(voff1),norm(voff2)))*2.0*expansion,3);
        srcface=[1 2 3];
    end
elseif(strcmp(cfg.srctype,'disk'))
    srcnode=orthdisk(cfg.srcpos,cfg.srcpos+cfg.srcdir,(cfg.srcparam1(1)*2.0)*expansion,3);
    srcface=delaunay(srcnode(:,1),srcnode(:,2));
else
    error(['source type not supported: ', cfg.srctype]);
end
