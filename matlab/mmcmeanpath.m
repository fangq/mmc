function avgpath=mmcmeanpath(detp,prop)
%
% avgpath=mmcmeanpath(detp,prop)
%
% Calculate the average pathlengths for each tissue type for a given source-detector pair
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     detp: the 2nd output from mmclab. detp can be either a struct or an array (detp.data)
%     prop: optical property list, as defined in the cfg.prop field of mmclab's input
%
% output:
%     avepath: the average pathlength for each tissue type 
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

detw=mmcdetweight(detp,prop);
avgpath=sum(detp.ppath.*repmat(detw(:),1,size(detp.ppath,2))) / sum(detw(:));
