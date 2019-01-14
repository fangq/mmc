function phi=sphdiffincident(r,theta,phi,cfg)
%
% phi=sphdiffincident(r,theta,phi,cfg)
%
% insident field of a sphere with a diffusion model
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     r,theta,phi: source position in spherical coordinates.
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
%     phi_inc=sphdiffincident(30,pi,0,cfg);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

% matlab's theta and phi are defined differently
[xs,ys,zs] = sph2cart(cfg.src(3),pi/2-cfg.src(2),cfg.src(1));
[x,y,z] = sph2cart(phi,pi/2-theta,r);
dist=sqrt((x-xs).*(x-xs)+(y-ys).*(y-ys)+(z-zs).*(z-zs));
phi=cfg.v./(4*pi*cfg.Dout*dist).*exp(i*cfg.kout*dist);

% if(isfield(cfg,'src2'))
%     [xs,ys,zs] = sph2cart(cfg.src2(3),pi/2-cfg.src2(2),cfg.src2(1));
%     dist=sqrt((x-xs).*(x-xs)+(y-ys).*(y-ys)+(z-zs).*(z-zs));
%     phi=phi-cfg.v./(4*pi*cfg.Dout*dist).*exp(i*cfg.kout*dist);
% end
