function res=sphdiffinterior(r,theta,phi,cfg)
res=0;
for l=0:cfg.maxl
  for m=-l:l
     res=res+(sphdiffCcoeff(m,l,cfg).*spbesselj(l,cfg.kin*r)).*spharmonic(l,m,theta,phi);
  end
end
