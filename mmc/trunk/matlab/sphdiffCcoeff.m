function C=sphdiffCcoeff(m,l,cfg)
%
% C=sphdiffCcoeff(m,l,cfg)
%
% sphere diffusion interior solution C coefficients
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     m: angular index
%     l: order
%     cfg: the problem domain setup, see sphdiffusioninfinite.m
%
% output:
%     res:  the coefficient at the specified order
%
% example:
%     C=sphdiffCcoeff(0,1,cfg)
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if((cfg.src(2)==pi | cfg.src(2)==0) & m~=0)
    C=0;
    return;
end
x=cfg.kout*cfg.a;
y=cfg.kin*cfg.a;
Dout=cfg.Dout;
Din=cfg.Din;
C=-i*cfg.v*cfg.kout/Dout*spbesselh(l,1,cfg.kout*cfg.src(1))*conj(spharmonic(l,m,cfg.src(2),cfg.src(3)))...
  *(Dout*x*spbesselh(l,1,x)*spbesseljprime(l,x)-Dout*x*spbesselhprime(l,1,x)*spbesselj(l,x))...
  /(Dout*x*spbesselhprime(l,1,x)*spbesselj(l,y)-Din*y*spbesselh(l,1,x)*spbesseljprime(l,y));
