function [Ja, newcfg]=mmcjmua(cfg,detp,seeds,detnum)
%
% Ja=mmcjmua(cfg,detp,seeds,detnum) (element based)
%
% Generate a time-domain Jacobian (sensitivity matrix) for absorption (mua) perturbation of a specified detector
% with configurations, detected photon seeds and detector readings from an initial MC simulation
%
% author: Ruoyang Yao (yaor <at> rpi.edu)
%         Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     cfg: the simulation configuration structure used for the initial MC simulation by mmclab
%     detp: detector readings from the initial MC simulation, must be a
%	        structure (supported after MMC v2016.4)
%     seeds: detected photon seeds from the initial MC simulation
%     detnum: the detector number whose detected photons will be replayed
%
% output:
%     Ja: a Jacobian for absorption perturbation of the specified detector
%		  number of rows is the number of the mesh nodes or elements
%		  number of columns is the number of time gates
%
% example:
%	  [cube,detp,ncfg,seeds]=mmclab(cfg);	% initial MC simulation
%	  Jmua = mmcjmua(ncfg,detp,seeds,1);    % generate absorption Jacobian of the first detector
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

% preprocessing of the mesh to get possible missing fields
newcfg=mmclab(cfg,'prep');
%newcfg.basisorder=0;

% specify the detector number to replay
newcfg.replaydet=detnum;

% select the photon seeds and data of the specified detector
newcfg.seed=seeds.data(:,detp.detid==newcfg.replaydet);
detp.data=detp.data(:,detp.detid==newcfg.replaydet);

% calculate the detected photon weight and arrival time
newcfg.replayweight=mmcdetweight(detp.data,newcfg.prop);
newcfg.replaytime=mmcdettime(detp.data,newcfg.prop);

% specify output type
newcfg.isnormalized=1; % should we leave this for users to decide?
newcfg.outputtype='wl';

% now replay detected photons
[jacob,detp1,newcfg]=mmclab(newcfg);

% validate if the replay is successful
if(all(ismember(round(detp1.data'*1e10)*1e-10,round(detp.data'*1e10)*1e-10,'rows')))
	% disp('replay is successful :-)');
	Ja=-jacob.data;
else
	error('replay failed :-(');
end
