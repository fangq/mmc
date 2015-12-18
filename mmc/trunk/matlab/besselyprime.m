function yp=besselyprime(n,z)
%
% yp=besselyprime(n,z)
%
% Neumann function (Bessel second kind) first order derivative 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Hankel function
%     z: input variable
%
% output:
%     yp: Neumann function (Bessel second kind) first order derivative 
%
% example:
%     yp=besselyprime(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%


yp=bessely(n-1,z)-n/z.*bessely(n,z);
