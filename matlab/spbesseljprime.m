function jp=spbesseljprime(n,z)
%
% jp=spbesseljprime(n,z)
%
% spherical Bessel function first order derivative
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Bessel function
%     z: input variable
%
% output:
%     jp: spherical Bessel function first order derivative
%
% example:
%     jp=spbesseljprime(0,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

jp=besseljprime(n+1/2,z).*sqrt(pi/(2*z))-sqrt(pi/2)*besselj(n+1/2,z)./(2*z.*sqrt(z));
