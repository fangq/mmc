function [Jmua, Jmus]=mmcjacobian(cfg,detp,seeds,detnum)
%
% [Jmua, Jmus]=mmcjacobian(cfg,detp,seeds,detnum) (element based)
%
% Generate time-domain Jacobians (sensitivity matrix) for absorption and scattering (mus) 
% perturbation of a specified source-detector pair with configurations, detected photon 
% seeds and detector readings from an initial MC simulation
%
% author: Ruoyang Yao (yaor <at> rpi.edu)
%         Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation configuration structure used for the initial MC simulation by mmclab
%     detp: detector readings from the initial MC simulation
%     seeds: detected photon seeds from the initial MC simulation
%     detnum: the detector number whose detected photons will be replayed
%
% output:
%     Jmua: the Jacobian for absorption coefficient for a specified source-detector pair
%           also output Jmua as a by-product.
%     Jmus: (optional) the Jacobian for scattering perturbation of a specified source detector pair
%		  number of rows is the number of the mesh elements
%		  number of columns is the number of time gates
%
% example:
%	  [cube,detp,ncfg,seeds]=mmclab(cfg);            % initial MC simulation
%	  [Jmua,Jmus] = mmcjacobian(ncfg,detp,seeds,1);  % generate scattering Jacobian of the first detector
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if(nargout==1)
    Jmua=mmcjmua(cfg,detp,seeds,detnum);
elseif(nargout==2)
    [Jmus,Jmua]=mmcjmus(cfg,detp,seeds,detnum);
end
