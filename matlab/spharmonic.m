function Y=spharmonic(l,m,theta,phi)
%
% Y=spharmonic(l,m,theta,phi)
%
% spherical harmonic function Y_{m,l}(theta,phi)
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     m: angular index
%     l: order
%     theta,phi: spherical angular coordinates, can be vectors
%
% output:
%     Y:  the values for Y_{m,l}(theta,phi)
%
% example:
%     Y=spharmonic(l,m,theta,phi)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

coeff=1;
oldm=m;
if(m<0)
    coeff=(-1)^m*prod(1:(l-m))/prod(1:(l+m));
    m=-m;
end
Lmn=legendre(l,cos((theta(:))'));
if l~=0
  Lmn=squeeze(Lmn(m+1,:))';
end
Lmn=reshape(Lmn,size(theta));
m=oldm;
Y=coeff*sqrt((2*l+1)*prod(1:(l-m))/(prod(1:(l+m))*4*pi)) ...
  *Lmn.*exp(i*m*phi);
