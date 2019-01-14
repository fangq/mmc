function [Jmus, Jmua]=mmcjmus(cfg,detp,seeds,detnum)
%
% [Jmus, Jmua]=mmcjmus(cfg,detp,seeds,detnum) (element based)
%
% Generate a time-domain Jacobian (sensitivity matrix) for scattering (mus) perturbation of a specified detector
% with configurations, detected photon seeds and detector readings from an initial MC simulation
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
%     Jmus: the Jacobian for scattering perturbation of a specified source detector pair
%		  number of rows is the number of the mesh elements
%		  number of columns is the number of time gates
%     Jmua: (optional) because calculating Jmus requires the Jacobian for mua, so one can
%           also output Jmua as a by-product.
%
% example:
%	  [cube,detp,ncfg,seeds]=mmclab(cfg);   % initial MC simulation
%	  Jmus = mmcjmus(ncfg,detp,seeds,1);    % generate scattering Jacobian of the first detector
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

[Jmua, newcfg]=mmcjmua(cfg,detp,seeds,detnum);

% specify output type 2
newcfg.outputtype='wp';

% replay detected photons for weighted partialpath
[jacob,detp2]=mmclab(newcfg);

% generate a map for scattering coefficient
elemp=newcfg.elemprop;
idx=find(elemp>0);
elemp(elemp<0)=0;
musmap0=newcfg.prop(elemp+1,2);

if(newcfg.basisorder==1)
    dim=4;
    nodemus=zeros(size(newcfg.node,1),1);
    nodevol=zeros(size(newcfg.node,1),1);
    for i=1:length(idx)
       nodemus(newcfg.elem(idx(i),1:dim))=nodemus(newcfg.elem(idx(i),1:dim))+newcfg.evol(idx(i))*musmap0(idx(i));
       nodevol(newcfg.elem(idx(i),1:dim))=nodevol(newcfg.elem(idx(i),1:dim))+newcfg.evol(idx(i));
    end
    musmap=nodemus./nodevol;
    musmap(isnan(musmap))=newcfg.prop(1,2);
else
    musmap=musmap0;
end

% divide by local scattering coefficient
jacob.data=jacob.data./repmat(musmap(:),1,size(jacob.data,2));
jacob.data(musmap<0)=0;

% validate if the replay is successful

if(all(ismember(round(detp2.data'*1e10)*1e-10,round(detp.data'*1e10)*1e-10,'rows')))
	%disp('replay is successful :-)');
	Jmus=Jmua+jacob.data;
else
	error('replay failed :-(');
end
