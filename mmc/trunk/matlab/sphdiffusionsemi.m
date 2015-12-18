function [res,xi,yi,zi] = sphdiffusionsemi(Reff,xrange,yrange,zrange,cfg)
%
% [res,xi,yi,zi]= sphdiffusionsemi(Reff,xrange,yrange,zrange,cfg)
%
% diffusion solution for a sphere inside a semi-infinite homogeneous medium 
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     Reff:  the effective reflection coeff.
%     xrange,yrange,zrange: a vector from where a grid will be created
%       		    and the phi values will be calculated
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
%   [phi,xi,yi,zi]=sphdiffusionsemi(0,-30:0.8:30,0,-30:0.8:30);
%   contourf(xi,zi,log10(abs(phi)),40)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%


if(nargin<5)
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
D = 1/(3*(cfg.omua+cfg.omusp));

zb = (1+Reff)/(1-Reff)*2*D;

z0 = 1/(cfg.omusp+cfg.omua);

src0=cfg.src;
cfg.src(1)=src0(1)-z0;

% real source full field for the real sphere
[res,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg); % S1,O1

cfg.src=src0;
cfg.src(2)=pi;
cfg.src(1)=src0(1)+z0+2*zb;
% image source full field for the real sphere
res2=sphdiffusioninfinite(xrange,yrange,zrange,cfg); % S2,O1

res=res-res2;

[P,T,R]=cart2sph(xi(:),yi(:),zi(:)); 
T=pi/2-T;  % matlab's theta and phi are defined differently
idx=find(R>cfg.a);

zrange=zrange+2*(src0(1)+zb);

cfg.src=src0;
cfg.src(2)=0;
cfg.src(1)=src0(1)+z0+2*zb;
% real source scattered field for the imaged sphere outside the real sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange,cfg); % S1,O2

res(idx)=res(idx)+res2(idx);

cfg.src=src0;
cfg.src(2)=0;
cfg.src(1)=src0(1)-z0;
% image source scattered field for the imaged sphere outside the real
% sphere
res2=sphdiffusionscatteronly(xrange,yrange,zrange,cfg); % S2,O2

res(idx)=res(idx)-res2(idx);

% high order terms are ignored, for example, the Phi(S2,O2) scattered by
% real sphere O1, or Phi(S1,O1) scatted by the mirrored sphere O2 etc
%
% 1st order approximation works only when the sphere is far away from
% source and interface so that multiple-scattering is weak.

%                      ____
%                    .'    `. real sphere
%                   /   O1   \
%                   |    o----|-----> x
%                   \    | a /  phi_interior=phi_interior(S1,O1)
%                    `.__|_.'               +phi_interior(S2,O1)
%                        |  phi_ext=phi_incident(S1,O1)+phi_scatter(S1,O1)
%                        |     - phi_incident(S2,O2)-phi_scatter(S2,O2)
%                        |     + phi_scatter(S1,O2) -phi_scatter(S2,O2)
%                        |     + high order terms(phi_scatter(S1,O1,O2),..)
%                     S1 x  real source
% _______________________|____________________z0______ true boundary
%                        |                    zb
% -----------------------|---------------------- extrapolated boundary
%           ^            |                    zb
%           |            |                    z0
%           |         S2 x  mirrored source  ---  S2=-S1
%           |            |
%           | src0(1)    |
%           |            |
%           |          ..|.  mirrored sphere
%           |        .'  | `.
%           v       /    | a \  (origin for the mirrored field)
%           ------- :    o----:-----> x
%                   \   O2   /
%                    `......'
%
%

