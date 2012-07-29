function [flux,detphoton]=mmclab(cfg)
%
%====================================================================
%      MMCLAB - Mesh-based Monte Carlo (MMC) for MATLAB/GNU Octave
%--------------------------------------------------------------------
%Copyright (c) 2012 Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%                   URL: http://mcx.sf.net/mmc/
%====================================================================
%
% Format:
%    [flux,detphoton]=mmclab(cfg);
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
%      cfg.debuglevel:  debug flag string, a subset of [MCBWDIOXATRPE]
%      cfg.nout:        [1.0] refractive index for medium type 0 (background)
%      cfg.minenergy:   terminate photon when weight less than this level (float) [0.0]
%      cfg.roulettesize:[10] size of Russian roulette
%      cfg.e0:          element that encloses the source (calc. if missing)
%      cfg.srctype:     source type, can be ["pencil"],"isotropic" or "cone"
%      cfg.srcparam:    1x4 vector for additional source parameter
%      cfg.unitinmm:    defines the unit in the input mesh [1.0]
%      
%
%      fields with * are required; options in [] are the default values
%
% Output:
%      flux: a struct array, with a length equals to that of cfg.
%            For each element of flux, flux(i).data is a 1D vector with
%            dimensions [size(cfg.node,1) total-time-gates] if cfg.basisorder=1,
%            or [size(cfg.elem,1) total-time-gates] if cfg.basisorder=0. 
%            The content of the array is the normalized flux (or others 
%            depending on cfg.outputtype) at each mesh node and time-gate.
%      detphoton: a struct array, with a length equals to that of cfg.
%            For each element of detphoton, detphoton(i).data is a 2D array with
%            dimensions [size(cfg.prop,1)+1 saved-photon-num]. The first row
%            is the ID(>0) of the detector that captures the photon; the second row
%            is the number of scattering events of the exitting photon; the rest rows
%            are the partial path lengths (in grid unit) traveling in medium 1 up 
%            to the last. If you set cfg.unitinmm, you need to multiply the path-lengths
%            to convert them to mm unit.
%
% Example:
%      cfg.nphoton=1e6;
%      [cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
%      cfg.elemprop=ones(size(cfg.elem,1),1);
%      cfg.srcpos=[30 30 0];
%      cfg.srcdir=[0 0 1];
%      cfg.prop=[0 0 1 1;0.005 1 0 1.37];
%      cfg.tstart=0;
%      cfg.tend=5e-9;
%      cfg.tstep=5e-10;
%      cfg.debuglevel='P';
%      % calculate the flux distribution with the given config
%      flux=mmclab(cfg);
%
%      cfgs(1)=cfg;
%      cfgs(2)=cfg;
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

len=length(cfg);
for i=1:len
    if(~isfield(cfg(i),'node') || ~isfield(cfg(i),'elem'))
        error('cfg.node or cfg.elem is missing');
    end
    if(~isfield(cfg(i),'elemprop') && size(cfg(i).elem,2)>4)
        cfg(i).elemprop=cfg(i).elem(:,5);
    end
    cfg(i).elem=meshreorient(cfg(i).node,cfg(i).elem(:,1:4));
    if(~isfield(cfg(i),'facenb'))
        cfg(i).facenb=faceneighbors(cfg(i).elem);
    end
    if(~isfield(cfg(i),'evol'))
        cfg(i).evol=elemvolume(cfg(i).node,cfg(i).elem);
    end
    if(~isfield(cfg(i),'srcpos'))
        error('cfg.srcpos field is missing');
    end
    if(~isfield(cfg(i),'srcdir'))
        error('cfg.srcdir field is missing');
    end
    if(~isfield(cfg(i),'e0'))
        cfg(i).e0=tsearchn(cfg(i).node,cfg(i).elem,cfg(i).srcpos);
    end
    if(isnan(cfg(i).e0))
        error('source is not located inside the mesh');
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

[flux,detphoton]=mmc(cfg);
