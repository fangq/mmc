function B=sphdiffBcoeff(m,l,cfg)
%
% B=sphdiffBcoeff(m,l,cfg)
%
% sphere diffusion exterior solution B coefficients
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     m: angular index
%     l: order
%     cfg: the problem domain setup, see sphdiffusioninfinite.m
%
% output:
%     res:  the coefficient at the specified order
%
% example:
%     B=sphdiffBcoeff(0,1,cfg)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

B=i*sphdiffAcoeff(m,l,cfg);
