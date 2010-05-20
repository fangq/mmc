function res=sphdiffexterior(r,theta,phi,cfg)
res=sphdiffincident(r,theta,phi,cfg)+sphdiffscatter(r,theta,phi,cfg);