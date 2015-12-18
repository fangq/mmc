function [T,P,R]=cart2sphorigin(xi,yi,zi,x0,y0,z0)
%
% [T,P,R]=cart2sphorigin(xi,yi,zi,x0,y0,z0)
%
% Converting Cartesian coordinates to spherical with a reset of origin
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     xi,yi,zi: list of Cartesian coordinates
%     x0,y0,z0: new origin in Cartesian coordinates
%
% output:
%     T: theta of the output spherical coordinates
%     P: phi of the output spherical coordinates
%     R: R of the output spherical coordinates
%
% example:
%     [T,P,R]=cart2sphorigin(1:5,1:5,1:5,2,2,2);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

xi=xi-x0;
yi=yi-y0;
zi=zi-z0;
[T,P,R]=cart2sph(xi(:),yi(:),zi(:));
