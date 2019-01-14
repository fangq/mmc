function [res,xi,yi,zi]=sphdiffusion(xrange,yrange,zrange,cfg)
%
% [res,xi,yi,zi]=sphdiffusion(xrange,yrange,zrange,cfg)
%
% diffusion solution for a sphere inside the infinite homogeneous medium 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     xrange,yrange,zrange: a vector from where a grid will be created
%       		    and the phi values will be calculated
%     cfg: the problem domain setup: 
%          cfg.v: speed of light in vacuum (mm/s)
%          cfg.a: sphere radius (mm)
%          cfg.omua: background (outside) mua (1/mm)
%          cfg.omusp: background (outside) mus' (1/mm)
%          cfg.imua: sphere (inside) mua (1/mm)
%          cfg.imusp: sphere (inside) mus' (1/mm)
%          cfg.src: spherical source position (R,theta,phi) R in mm
%          cfg.maxl: maximum serial expansion terms
%          cfg.omega: DPDW modulation frequency
%
% output:
%     res:  the output fluence for both interior and exterior regions
%
% example:
%   [phi_ana,xa,ya,za]=sphdiffusion(-30:0.8:30,0,-30:0.8:30);
%   contourf(xa,za,log10(abs(phi_ana)),40)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%


if(nargin<4)
	cfg.v=299792458000;
	cfg.a=10;
	cfg.omua=0.002;
	cfg.omusp=0.990;
	cfg.imua=0.050;
	cfg.imusp=0.500;
	cfg.src=[30,pi,0];
	cfg.maxl=20;
	cfg.omega=0;
end

cfg.Din=cfg.v/(3*cfg.imusp);
cfg.Dout=cfg.v/(3*cfg.omusp);
cfg.kin=sqrt((-cfg.v*cfg.imua+i*cfg.omega)/cfg.Din);
cfg.kout=sqrt((-cfg.v*cfg.omua+i*cfg.omega)/cfg.Dout);

[xi,yi,zi]=meshgrid(xrange,yrange,zrange);

[P,T,R]=cart2sph(xi(:),yi(:),zi(:)); % matlab's theta and phi are defined differently
T=pi/2-T;

idx=find(R>cfg.a);
res=zeros(length(R),1);
res(idx)=sphdiffexterior(R(idx),T(idx),P(idx),cfg);

idx=find(R<=cfg.a);
res(idx)=sphdiffinterior(R(idx),T(idx),P(idx),cfg);

res=squeeze(reshape(res,size(xi)));
xi=squeeze(xi);
yi=squeeze(yi);
zi=squeeze(zi);

%%-------------------------------------------------------

function res=sphdiffexterior(r,theta,phi,cfg)
res=sphdiffincident(r,theta,phi,cfg)+sphdiffscatter(r,theta,phi,cfg);

%%-------------------------------------------------------

function res=sphdiffinterior(r,theta,phi,cfg)
res=0;
for l=0:cfg.maxl
  for m=-l:l
     res=res+(sphdiffCcoeff(m,l,cfg).*spbesselj(l,cfg.kin*r)).*spharmonic(l,m,theta,phi);
  end
end

%%-------------------------------------------------------

function phi=sphdiffincident(r,theta,phi,cfg)
% matlab's theta and phi are defined differently
[xs,ys,zs] = sph2cart(cfg.src(3),pi/2-cfg.src(2),cfg.src(1));
[x,y,z] = sph2cart(phi,pi/2-theta,r);
dist=sqrt((x-xs).*(x-xs)+(y-ys).*(y-ys)+(z-zs).*(z-zs));
phi=cfg.v./(4*pi*cfg.Dout*dist).*exp(i*cfg.kout*dist);

% if(isfield(cfg,'src2'))
%     [xs,ys,zs] = sph2cart(cfg.src2(3),pi/2-cfg.src2(2),cfg.src2(1));
%     [x,y,z] = sph2cart(phi,pi/2-theta,r);
%     dist=sqrt((x-xs).*(x-xs)+(y-ys).*(y-ys)+(z-zs).*(z-zs));
%     phi=phi-cfg.v./(4*pi*cfg.Dout*dist).*exp(i*cfg.kout*dist);
% end

%%-------------------------------------------------------

function res=sphdiffscatter(r,theta,phi,cfg)
res=zeros(size(r));
for l=0:cfg.maxl
  for m=-l:l
     res=res+(sphdiffAcoeff(m,l,cfg)*spbesselj(l,cfg.kout*r)+sphdiffBcoeff(m,l,cfg)...
        *spbessely(l,cfg.kout*r)).*spharmonic(l,m,theta,phi);
  end
end

%%-------------------------------------------------------

function A=sphdiffAcoeff(m,l,cfg)
x=cfg.kout*cfg.a;
y=cfg.kin*cfg.a;
Dout=cfg.Dout;
Din=cfg.Din;
A=-i*cfg.v*cfg.kout/Dout*spbesselh(l,1,cfg.kout*cfg.src(1))*conj(spharmonic(l,m,pi,0))...
  *(Dout*x*spbesseljprime(l,x)*spbesselj(l,y)-Din*y*spbesselj(l,x)*spbesseljprime(l,y))...
  /(Dout*x*spbesselhprime(l,1,x)*spbesselj(l,y)-Din*y*spbesselh(l,1,x)*spbesseljprime(l,y));

%%-------------------------------------------------------

function B=sphdiffBcoeff(m,l,cfg)
B=i*sphdiffAcoeff(m,l,cfg);

%%-------------------------------------------------------

function C=sphdiffCcoeff(m,l,cfg)
x=cfg.kout*cfg.a;
y=cfg.kin*cfg.a;
Dout=cfg.Dout;
Din=cfg.Din;
C=-i*cfg.v*cfg.kout/Dout*spbesselh(l,1,cfg.kout*cfg.src(1))*conj(spharmonic(l,m,pi,0))...
  *(Dout*x*spbesselh(l,1,x)*spbesseljprime(l,x)-Dout*x*spbesselhprime(l,1,x)*spbesselj(l,x))...
  /(Dout*x*spbesselhprime(l,1,x)*spbesselj(l,y)-Din*y*spbesselh(l,1,x)*spbesseljprime(l,y));

%%-------------------------------------------------------

function Y=spharmonic(l,m,theta,phi)
coeff=1;
oldm=m;
if(m<0)
    coeff=(-1)^m*prod(1:(l-m))/prod(1:(l+m));
    m=-m;
end
Lmn=legendre(l,cos(theta));
if l~=0
  Lmn=squeeze(Lmn(m+1,:))';
end
m=oldm;
Y=coeff*sqrt((2*l+1)*prod(1:(l-m))/(prod(1:(l+m))*4*pi)) ...
  *Lmn.*exp(i*m*phi);

%%-------------------------------------------------------

function jn=spbesselj(n,z)
jn=besselj(n+1/2,z).*sqrt(pi./(2*z));

%%-------------------------------------------------------

function yn=spbessely(n,z)
yn=bessely(n+1/2,z).*sqrt(pi./(2*z));

%%-------------------------------------------------------

function hn=spbesselh(n,k,z)
if(k==1)
    hn=spbesselj(n,z)+i*spbessely(n,z);
elseif(k==2)
    hn=spbesselj(n,z)-i*spbessely(n,z);
else
    error('wrong value for the second parameter');
end
%%-------------------------------------------------------

function jp=besseljprime(n,z)
jp=besselj(n-1,z)-n/z.*besselj(n,z);

%%-------------------------------------------------------

function yp=besselyprime(n,z)
yp=bessely(n-1,z)-n/z.*bessely(n,z);

%%-------------------------------------------------------

function hp=besselhprime(n,k,z)
hp=besselh(n-1,k,z)-n/z.*besselh(n,k,z);

%%-------------------------------------------------------

function jp=spbesseljprime(n,z)
jp=besseljprime(n+1/2,z).*sqrt(pi/(2*z))-sqrt(pi/2)*besselj(n+1/2,z)./(2*z.*sqrt(z));

%%-------------------------------------------------------

function yp=spbesselyprime(n,z)
yp=besselyprime(n+1/2,z).*sqrt(pi/(2*z))-sqrt(pi/2)*bessely(n+1/2,z)./(2*z.*sqrt(z));

%%-------------------------------------------------------

function hn=spbesselhprime(n,k,z)
if(k==1)
    hn=spbesseljprime(n,z)+i*spbesselyprime(n,z);
elseif(k==2)
    hn=spbesseljprime(n,z)-i*spbesselyprime(n,z);
else
    error('wrong value for the second parameter');
end



