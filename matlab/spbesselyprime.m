function yp=spbesselyprime(n,z)
%
% yp=spbesselyprime(n,z)
%
% spherical Neumann function first order derivative 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Neumann function
%     z: input variable
%
% output:
%     yp:  spherical Neumann function first order derivative 
%
% example:
%     yp=spbesselyprime(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

yp=besselyprime(n+1/2,z).*sqrt(pi/(2*z))-sqrt(pi/2)*bessely(n+1/2,z)./(2*z.*sqrt(z));
