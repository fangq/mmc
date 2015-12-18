function yn=spbessely(n,z)
%
% yn=spbessely(n,z)
%
% spherical Neumann function 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Neumann function
%     z: input variable
%
% output:
%     yn:  spherical Neumann function first order derivative 
%
% example:
%     yn=spbessely(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

yn=bessely(n+1/2,z).*sqrt(pi./(2*z));
