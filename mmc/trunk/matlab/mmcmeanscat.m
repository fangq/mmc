function avgnscat=mmcmeanscat(detp,prop)
%
% avgnscat=mmcmeanscat(detp,prop)
%
% Calculate the average scattering event counts for each tissue type for a given source-detector pair
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     detp: the 2nd output from mmclab. detp can be either a struct or an array (detp.data)
%     prop: optical property list, as defined in the cfg.prop field of mmclab's input
%
% output:
%     avgnscat: the average scattering event count for each tissue type 
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

detw=mmcdetweight(detp,prop);
avgnscat=sum(detp.nscat.*repmat(detw(:),1,size(detp.ppath,2))) / sum(detw(:));
