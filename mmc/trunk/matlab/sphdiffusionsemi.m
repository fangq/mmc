function [res,xi,yi,zi] = sphdiffusionsemi(Reff,xrange,yrange,zrange,cfg)
%  Phi = sphdiffusionsemi(Reff,xrange,yrange,zrange,cfg)
%
%  semi-infinite medium analytical solution to diffusion model
%
%    author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
%    input:
%        mua:   the absorption coefficients in 1/mm
%        musp:  the reduced scattering coefficients in 1/mm
%        Reff:  the effective reflection coeff.
%        srcpos:array for the source positions (x,y,z)
%        detpos:array for the detector positions (x,y,z)
%
%    output:
%        Phi:  the output fluence for all source/detector pairs
%
%    this file is part of Mesh-based Monte Carlo (MMC)
%    License: GPLv3, see http://mcx.sf.net/?MMC for details


if(nargin<5)
	cfg.v=299792458000;
	cfg.a=20;
	cfg.omua=0.002;
	cfg.omusp=0.990;
	cfg.imua=0.050;
	cfg.imusp=0.500;
% 	cfg.imua=0.002;
% 	cfg.imusp=0.990;
	cfg.src=[30,pi,0];
	cfg.maxl=20;
	cfg.omega=0;
end

D = 1/(3*(cfg.omua+cfg.omusp));
zb = (1+Reff)/(1-Reff)*2*D;

z0 = 1/(cfg.omusp+cfg.omua);

src0=cfg.src;
cfg.src(1)=cfg.src(1)-z0;

[res,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg);

cfg.src=src0;
cfg.src(1)=cfg.src(1)+z0+2*zb;

[res2,xi,yi,zi]=sphdiffusioninfinite(xrange,yrange,zrange,cfg);

res=res-res2;