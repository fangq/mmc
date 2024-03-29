function res = sphdiffscatter(r, theta, phi, cfg)
%
% res=sphdiffscatter(r,theta,phi,cfg)
%
% sphere exterior scattered field
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
%     phi_scat=sphdiffscatter(30,pi,0,cfg);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

% if(cfg.src(2)<pi/2)
%        theta=pi-theta; % mirror up-to-down
% end
res = zeros(size(r));
for l = 0:cfg.maxl
    for m = -l:l
        res = res + (sphdiffAcoeff(m, l, cfg) * spbesselj(l, cfg.kout * r) + sphdiffBcoeff(m, l, cfg) * ...
                     spbessely(l, cfg.kout * r)) .* spharmonic(l, m, theta, phi);
    end
end
