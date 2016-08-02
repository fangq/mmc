function Ja=mmcmuaJacobian(cfg,detp,seeds,detnum)
%
% Ja=mmcmuaJacobian(cfg,detp,seeds,detnum) (element based)
%
% Generate a time-domain Jacobian (sensitivity matrix) for absorption (mua) perturbation of a specified detector
% with configurations, detected photon seeds and detector readings from an initial MC simulation
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%		  Ruoyang Yao (yaor <at> rpi.edu)
%
% input:
%     cfg: the simulation configuration structure used for the initial MC simulation by mmclab
%	  detp: detector readings from the initial MC simulation
%     seeds: detected photon seeds from the initial MC simulation
%     detnum: the detector number whose detected photons will be replayed
%
% output:
%     Ja: a Jacobian for absorption perturbation of the specified detector
%		  number of rows is the number of the mesh nodes or elements
%		  number of columns is the number of time gates
%
% example:
%	  [cube,detp,~,seeds]=mmclab(newcfg);		% initial MC simulation
%	  J_mua = mmcmuaJacobian(cfg,detp,seeds,1);	% generate absorption Jacobian of the first detector
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

% preprocessing of the mesh to get possible missing fields
newcfg=mmclab(cfg,'prep');
newcfg.basisorder=0;

% specify the detector number to replay
newcfg.replaydet=detnum;

% select the photon seeds and data of the specified detector
newcfg.seed=seeds.data(:,detp.data(1,:)==newcfg.replaydet);
detp.data=detp.data(:,detp.data(1,:)==newcfg.replaydet);

% calculate the detected photon weight and arrival time
newcfg.replayweight=mmcdetweight(detp.data,newcfg.prop);
newcfg.replaytime=mmcdettime(detp.data,newcfg.prop);

% specify output type
newcfg.isnormalized=1;
newcfg.outputtype='wl';

% now replay detected photons
[cube,detp1,~,~]=mmclab(newcfg);

% validate if the replay is successful
if(all(ismember(round(detp.data'*1e10)*1e-10,round(detp1.data'*1e10)*1e-10,'rows')))
	disp('replay is successful :-)');
	Ja=-cube.data;
else
	error('replay failed :-(');
end
