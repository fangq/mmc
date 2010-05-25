function res=sphdiffscatter(r,theta,phi,cfg)
%if(cfg.src(2)<pi/2)
%        theta=pi-theta; % mirror up-to-down
%end
res=zeros(size(r));
for l=0:cfg.maxl
  for m=-l:l
     res=res+(sphdiffAcoeff(m,l,cfg)*spbesselj(l,cfg.kout*r)+sphdiffBcoeff(m,l,cfg)...
        *spbessely(l,cfg.kout*r)).*spharmonic(l,m,theta,phi);
  end
end
