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
%    cfg: a struct, or struct array. Each element of cfg defines 
%         the parameters associated with a simulation. 
%
%    It may contain the following fields:
%
%     *cfg.nphoton:    the total number of photons to be simulated (integer)
%     *cfg.prop:       an N by 4 array, each row specifies [mua, mus, g, n] in order.
%                      the first row corresponds to medium type 0 which is 
%                      typically [0 0 1 1]. The second row is type 1, and so on.
%     *cfg.node:       node array for the input tetrahedral mesh, 3 columns: (x,y,z)
%     *cfg.elem:       element array for the input tetrahedral mesh, 4 columns
%     *cfg.tstart:     starting time of the simulation (in seconds)
%     *cfg.tstep:      time-gate width of the simulation (in seconds)
%     *cfg.tend:       ending time of the simulation (in second)
%     *cfg.srcpos:     a 1 by 3 vector, the position of the source in grid unit
%     *cfg.srcdir:     a 1 by 3 vector, specifying the incident vector
%      cfg.session:    a string for output file names (used when no return variables)
%      cfg.seed:       seed for the random number generator (integer) [0]
%      cfg.detpos:     an N by 4 array, each row specifying a detector: [x,y,z,radius]
%      cfg.isreflect:  [1]-consider refractive index mismatch, 0-matched index
%      cfg.isrefint:   1-ref. index mismatch at inner boundaries, [0]-matched index
%      cfg.isnormalized:[1]-normalize the output flux to unitary source, 0-no reflection
%      cfg.minenergy:  terminate photon when weight less than this level (float) [0.0]
%      cfg.unitinmm:   defines the length unit for a grid edge length [1.0]
%
%      fields with * are required; options in [] are the default values
%
% Output:
%      flux: a struct array, with a length equals to that of cfg.
%            For each element of flux, flux(i).data is a 4D array with
%            dimensions specified by [size(vol) total-time-gates]. 
%            The content of the array is the normalized flux at 
%            each voxel of each time-gate.
%      detphoton: a struct array, with a length equals to that of cfg.
%            For each element of detphoton, detphoton(i).data is a 2D array with
%            dimensions [size(cfg.prop,1)+1 saved-photon-num]. The first row
%            is the ID(>0) of the detector that captures the photon; the second row
%	     is the number of scattering events of the exitting photon; the rest rows
%	     are the partial path lengths (in grid unit) traveling in medium 1 up 
%            to the last. If you set cfg.unitinmm, you need to multiply the path-lengths
%            to convert them to mm unit.
%
%      if detphoton is ignored, the detected photon will be saved in a .mch file 
%      if cfg.issavedeet=1; if no output is given, the flux will be saved to a 
%      .mc2 file if cfg.issave2pt=1 (which is true by default).
%
% Example:
%      cfg.nphoton=1e7;
%      cfg.vol=uint8(ones(60,60,60));
%      cfg.srcpos=[30 30 1];
%      cfg.srcdir=[0 0 1];
%      cfg.prop=[0 0 1 1;0.005 1 0 1.37];
%      cfg.tstart=0;
%      cfg.tend=5e-9;
%      cfg.tstep=5e-10;
%      % calculate the flux distribution with the given config
%      flux=mmclab(cfg);
%
%      cfgs(1)=cfg;
%      cfgs(2)=cfg;
%      cfgs(1).isreflect=0;
%      cfgs(2).isreflect=1;
%      cfgs(2).issavedet=1;
%      cfgs(2).detpos=[30 20 1 1;30 40 1 1;20 30 1 1;40 30 1 1];
%      % calculate the flux and partial path lengths for the two configurations
%      [fluxs,detps]=mmclab(cfgs);
%
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
