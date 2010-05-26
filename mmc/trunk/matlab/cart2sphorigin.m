function [T,P,R]=cart2sphorigin(xi,yi,zi,x0,y0,z0)
xi=xi-x0;
yi=yi-y0;
zi=zi-z0;
[T,P,R]=cart2sph(xi(:),yi(:),zi(:));
