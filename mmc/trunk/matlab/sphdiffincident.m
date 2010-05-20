function phi=sphdiffincident(r,theta,phi,cfg)
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