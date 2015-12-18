function jn=spbesselj(n,z)
%
% jn=spbesselj(n,z)
%
% spherical Bessel function
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Bessel function
%     z: input variable
%
% output:
%     jn: spherical Bessel function first order derivative
%
% example:
%     jn=spbesselj(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

jn=besselj(n+1/2,z).*sqrt(pi./(2*z));
