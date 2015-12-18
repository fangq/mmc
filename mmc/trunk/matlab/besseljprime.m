function jp=besseljprime(n,z)
%
% jp=besseljprime(n,z)
%
% Bessel function (Bessel first kind) first order derivative 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Hankel function
%     z: input variable
%
% output:
%     jp: Bessel function (Bessel first kind) first order derivative
%
% example:
%     jp=besseljprime(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

jp=besselj(n-1,z)-n/z.*besselj(n,z);
