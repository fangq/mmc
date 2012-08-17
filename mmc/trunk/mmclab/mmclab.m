function varargout=mmclab(cfg,type)
%
%====================================================================
%      MMCLAB - Mesh-based Monte Carlo (MMC) for MATLAB/GNU Octave
%--------------------------------------------------------------------
%Copyright (c) 2012 Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%                   URL: http://mcx.sf.net/mmc/
%====================================================================
%
% Format:
%    [flux,detphoton,ncfg]=mmclab(cfg,type);
%
% Input:
%    cfg: a struct, or struct array. Each element in cfg defines 
%         a set of parameters for a simulation. 
%
%    It may contain the following fields:
%
%     *cfg.nphoton:     the total number of photons to be simulated (integer)
%     *cfg.prop:        an N by 4 array, each row specifies [mua, mus, g, n] in order.
%                       the first row corresponds to medium type 0 which is 
%                       typically [0 0 1 1]. The second row is type 1, and so on.
%     *cfg.node:        node array for the input tetrahedral mesh, 3 columns: (x,y,z)
%     *cfg.elem:        element array for the input tetrahedral mesh, 4 columns
%     *cfg.elemprop:    element property index for input tetrahedral mesh
%     *cfg.tstart:      starting time of the simulation (in seconds)
%     *cfg.tstep:       time-gate width of the simulation (in seconds)
%     *cfg.tend:        ending time of the simulation (in second)
%     *cfg.srcpos:      a 1 by 3 vector, the position of the source in grid unit
%     *cfg.srcdir:      a 1 by 3 vector, specifying the incident vector
%     -cfg.facenb:      element face neighbohood list (calculated by faceneighbors())
%     -cfg.evol:        element volume (calculated by elemvolume() with iso2mesh)
%     -cfg.e0:          the element ID enclosing the source, if not defined,
%                       it will be calculated by tsearchn(node,elem,srcpos);
%                       if cfg.e0 is set as one of the following characters,
%                       mmclab will do an initial ray-tracing and move
%                       srcpos to the first intersection to the surface:
%                       '>': search along the forward (srcdir) direction
%                       '<': search along the backward direction
%                       '-': search both directions
%      cfg.seed:        seed for the random number generator (integer) [0]
%      cfg.detpos:      an N by 4 array, each row specifying a detector: [x,y,z,radius]
%      cfg.isreflect:   [1]-consider refractive index mismatch, 0-matched index
%      cfg.isnormalized:[1]-normalize the output flux to unitary source, 0-no reflection
%      cfg.isspecular:  [1]-calculate specular reflection if source is outside
%      cfg.basisorder:  [1]-linear basis, 0-piece-wise constant basis
%      cfg.outputformat:['ascii'] or 'bin' (in 'double')
%      cfg.outputtype:  [X] - output flux, F - fluence, E - energy deposit
%      cfg.method:      ray-tracing method, [P]:Plucker, H:Havel (SSE4),
%                       B: partial Badouel, S: branchless Badouel (SSE)
%      cfg.debuglevel:  debug flag string, a subset of [MCBWDIOXATRPE], no space
%      cfg.nout:        [1.0] refractive index for medium type 0 (background)
%      cfg.minenergy:   terminate photon when weight less than this level (float) [0.0]
%      cfg.roulettesize:[10] size of Russian roulette
%      cfg.srctype:     source type, can be ["pencil"],"isotropic" or "cone"
%      cfg.srcparam:    1x4 vector for additional source parameter
%      cfg.unitinmm:    defines the unit in the input mesh [1.0]
%
%      fields marked with * are required; options in [] are the default values
%      fields marked with - are calculated if not given (can be faster if precomputed)
%
%    type: omit or 'omp' for multi-threading version; 'sse' for the SSE4 MMC,
%          the SSE4 version is about 25% faster, but requires newer CPUs.
%
% Output:
%      flux: a struct array, with a length equals to that of cfg.
%            For each element of flux, flux(i).data is a 1D vector with
%            dimensions [size(cfg.node,1) total-time-gates] if cfg.basisorder=1,
%            or [size(cfg.elem,1) total-time-gates] if cfg.basisorder=0. 
%            The content of the array is the normalized flux (or others 
%            depending on cfg.outputtype) at each mesh node and time-gate.
%      detphoton: (optional) a struct array, with a length equals to that of cfg.
%            For each element of detphoton, detphoton(i).data is a 2D array with
%            dimensions [size(cfg.prop,1)+1 saved-photon-num]. The first row
%            is the ID(>0) of the detector that captures the photon; the second row
%            is the number of scattering events of the exitting photon; the rest rows
%            are the partial path lengths (in grid unit) traveling in medium 1 up 
%            to the last. If you set cfg.unitinmm, you need to multiply the path-lengths
%            to convert them to mm unit.
%      ncfg: (optional), if given, mmclab returns the preprocessed cfg structure,
%            including the calculated subfields (marked by "-"). This can be
%            used as the input to avoid repetitive preprocessing.
%
% Example:
%      cfg.nphoton=1e5;
%      [cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
%      cfg.elemprop=ones(size(cfg.elem,1),1);
%      cfg.srcpos=[30 30 0];
%      cfg.srcdir=[0 0 1];
%      cfg.prop=[0 0 1 1;0.005 1 0 1.37];
%      cfg.tstart=0;
%      cfg.tend=5e-9;
%      cfg.tstep=5e-10;
%      cfg.debuglevel='TP';
%      % calculate the flux distribution with the given config
%      [flux detp ncfg]=mmclab(cfg);
%
%      cfgs(1)=ncfg;
%      cfgs(2)=ncfg;
%      cfgs(1).isreflect=0;
%      cfgs(2).isreflect=1;
%      cfgs(2).detpos=[30 20 0 1;30 40 0 1;20 30 1 1;40 30 0 1];
%      % calculate the flux and partial path lengths for the two configurations
%      [fluxs,detps]=mmclab(cfgs);
%
%
% This function is part of Mesh-based Monte Carlo (MMC) URL: http://mcx.sf.net/mmc/
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

if(nargin==0)
    error('input field cfg must be defined');
end
if(~isstruct(cfg))
    error('cfg must be a struct or struct array');
end

len=length(cfg);
for i=1:len
    if(~isfield(cfg(i),'node') || ~isfield(cfg(i),'elem'))
        error('cfg.node or cfg.elem is missing');
    end
    if(~isfield(cfg(i),'elemprop') && size(cfg(i).elem,2)>4)
        cfg(i).elemprop=cfg(i).elem(:,5);
    end
    cfg(i).elem=meshreorient(cfg(i).node,cfg(i).elem(:,1:4));
    if(~isfield(cfg(i),'facenb') || isempty(cfg(i).facenb))
        cfg(i).facenb=faceneighbors(cfg(i).elem);
    end
    if(~isfield(cfg(i),'evol') || isempty(cfg(i).evol))
        cfg(i).evol=elemvolume(cfg(i).node,cfg(i).elem);
    end
    if(~isfield(cfg(i),'srcpos'))
        error('cfg.srcpos field is missing');
    end
    if(~isfield(cfg(i),'srcdir'))
        error('cfg.srcdir field is missing');
    end
    if(~isfield(cfg(i),'e0') || isempty(cfg(i).e0))
        cfg(i).e0=tsearchn(cfg(i).node,cfg(i).elem,cfg(i).srcpos);
    end
    if(isnan(cfg(i).e0) || ischar(cfg(i).e0))
        disp('searching initial element ...');
        face=volface(cfg(i).elem);
        [t,u,v,idx]=raytrace(cfg(i).srcpos,cfg(i).srcdir,cfg(i).node,face);
        if(isempty(idx))
            error('ray does not intersect with the mesh');
        else
            t=t(idx);
            if(cfg(i).e0=='>')
                idx1=find(t>=0);
            elseif(cfg(i).e0=='<')
                idx1=find(t<=0);
            elseif(isnan(cfg(i).e0) || cfg(i).e0=='-')
                idx1=1:length(t);
            else
                error('ray direction specifier is not recognized');
            end
            if(isempty(idx1))
                error('no intersection is found along the ray direction');
            end
            t0=abs(t(idx1));
            [tmin,loc]=min(t0);
            faceidx=idx(idx1(loc));

            % update source position
            cfg(i).srcpos=cfg(i).srcpos+t(idx1(loc))*cfg(i).srcdir;

            % find initial element id
            felem=sort(face(faceidx,:));
            f=cfg(i).elem;
            f=[f(:,[1,2,3]);
               f(:,[2,1,4]);
               f(:,[1,3,4]);
               f(:,[2,4,3])];
            [tf,loc]=ismember(felem,sort(f,2),'rows');
            loc=mod(loc,size(cfg(i).elem,1));
            if(loc==0) loc=size(cfg(i).elem,1); end
            cfg(i).e0=loc;
        end
    end
    if(~isfield(cfg(i),'elemprop'))
        error('cfg.elemprop field is missing');
    end
    if(~isfield(cfg(i),'nphoton'))
        error('cfg.nphoton field is missing');
    end
    if(~isfield(cfg(i),'prop') || size(cfg(i).prop,1)<max(cfg(i).elemprop)+1)
        error('cfg.prop field is missing or insufficient');
    end
end

% must do a fflush, otherwise octave buffers the output until complete
if(exist('OCTAVE_VERSION'))
   fflush(stdout);
end

mmcout=nargout;
if(nargout>=3)
    mmcout=2;
    varargout{3}=cfg;
end

if(nargin<2)
  [varargout{1:mmcout}]=mmc(cfg);
elseif(strcmp(type,'omp'))
  [varargout{1:mmcout}]=mmc(cfg);
elseif(strcmp(type,'sse'))
  [varargout{1:mmcout}]=mmc_sse(cfg);
else
  error('type is not recognized');
end
