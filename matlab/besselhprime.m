function hp=besselhprime(n,k,z)
%
% hp=besselhprime(n,k,z)
%
% Hankel function first order derivative 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Hankel function
%     k: kind of the Hankel function
%     z: input variable
%
% output:
%     hn: Hankel function first order derivative 
%
% example:
%     hn=besselhprime(0,1,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

hp=besselh(n-1,k,z)-n/z.*besselh(n,k,z);
