function [res,xi,yi,zi]=sphdiffusionscatteronly(xrange,yrange,zrange,cfg)

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
res(idx)=sphdiffscatter(R(idx),T(idx),P(idx),cfg);

res=squeeze(reshape(res,size(xi)));
xi=squeeze(xi);
yi=squeeze(yi);
zi=squeeze(zi);