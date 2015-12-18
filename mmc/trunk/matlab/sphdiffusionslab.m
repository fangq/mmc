function [res,xi,yi,zi] = sphdiffusionslab(Reff1,Reff2,h,xrange,yrange,zrange,cfg)
%
% [res,xi,yi,zi]= sphdiffusionslab(Reff,h,xrange,yrange,zrange,cfg)
%
% diffusion solution for a sphere inside an infinite homogeneous slab 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Reff:  the effective reflection coeff.
%     xrange,yrange,zrange: a vector from where a grid will be created
%       		    and the phi values will be calculated
%     h: the height of the slab
%     cfg: domain structure for internal/external parameters
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
%     res:  the output fluence for both the interior and exterior
%     regions
%
% example:
%   [phi_ana,xa,ya,za]=sphdiffusionslab(0,0,60,-30:0.8:30,0,-30:0.8:30);
%   contourf(xa,za,log10(abs(phi_ana)),40)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if(nargin<7)
	cfg.v=299792458000;  % mm/s
	cfg.a=10;            % radius, mm
	cfg.omua=0.002;      % outside mua 1/mm
	cfg.omusp=0.990;     % outside mus' 1/mm
	cfg.imua=0.050;
	cfg.imusp=0.500;
% 	cfg.imua=0.002;
% 	cfg.imusp=0.990;
	cfg.src=[30,pi,0];
	cfg.maxl=20;
	cfg.omega=0;
end
cfg.Din=cfg.v/(3*cfg.imusp);
cfg.Dout=cfg.v/(3*cfg.omusp);
cfg.kin=sqrt((-cfg.v*cfg.imua+i*cfg.omega)/cfg.Din);
cfg.kout=sqrt((-cfg.v*cfg.omua+i*cfg.omega)/cfg.Dout);
src0=cfg.src;

[res,xi,yi,zi]=sphdiffusionsemi(Reff1,xrange,yrange,zrange,cfg);

D = 1/(3*(cfg.omua+cfg.omusp));
zb = (1+Reff2)/(1-Reff2)*2*D;
z0 = 1/(cfg.omusp+cfg.omua);

cfg.src=[2*h-src0(1)+2*zb-z0, pi-src0(2), src0(3)];

% image source at the upper interface full field for the real sphere
[res2,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg); % S2,O1

res=res-res2;

cfg.src=[2*h-src0(1)+2*zb-z0+2*(z0+zb), pi-src0(2), src0(3)];

% image source at the upper interface full field for the real sphere
[res2,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg); % S2,O1

res=res+res2;


[P,T,R]=cart2sph(xi(:),yi(:),zi(:));
T=pi/2-T;  % matlab's theta and phi are defined differently
idx=find(R>cfg.a);


cfg.src=[2*h-src0(1)-z0, src0(2), src0(3)];
% real source scattered field for the imaged sphere outside the real sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange-2*(h-src0(1)+zb),cfg); % S1,O2

res(idx)=res(idx)+res2(idx);

cfg.src=[2*h-src0(1)+2*zb+z0, src0(2), src0(3)]; 
% image source scattered field for the imaged sphere outside the real
% sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange-2*(h-src0(1)+zb),cfg); % S2,O2

res(idx)=res(idx)-res2(idx);


% this translate the grid to the origin of the mirrored sphere

cfg.src=[src0(1)-z0, pi-src0(2), src0(3)];
% real source scattered field for the imaged sphere outside the real sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange-2*(h-src0(1)+zb),cfg); % S1,O2

res(idx)=res(idx)-res2(idx);

cfg.src=[src0(1)+2*zb+z0, pi-src0(2), src0(3)];
% image source scattered field for the imaged sphere outside the real
% sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange-2*(h-src0(1)+zb),cfg); % S2,O2

res(idx)=res(idx)+res2(idx);


% high order terms are ignored, for example, the Phi(S2,O2) scattered by
% real sphere O1, or Phi(S1,O1) scatted by the mirrored sphere O2 etc
%
% 1st order approximation works only when the sphere is far away from
% source and interface so that multiple-scattering is weak.
