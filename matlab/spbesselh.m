function hn=spbesselh(n,k,z)
%
% hn=spbesselhprime(n,k,z)
%
% spherical Hankel function 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     n: order of the spherical Hankel function
%     k: kind of the Hankel function
%     z: input variable
%
% output:
%     hn: spherical Hankel function 
%
% example:
%     hn=spbesselh(0,1,1)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if(k==1)
    hn=spbesselj(n,z)+i*spbessely(n,z);
elseif(k==2)
    hn=spbesselj(n,z)-i*spbessely(n,z);
else
    error('wrong value for the second parameter');
end
