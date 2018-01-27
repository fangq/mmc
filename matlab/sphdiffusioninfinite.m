function [res,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg)
%
% [res,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg)
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
%   [phi_ana,xa,ya,za]=sphdiffusioninfinite(-30:0.8:30,0,-30:0.8:30);
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
