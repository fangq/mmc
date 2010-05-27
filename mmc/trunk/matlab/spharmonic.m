function Y=spharmonic(l,m,theta,phi)
coeff=1;
oldm=m;
if(m<0)
    coeff=(-1)^m*prod(1:(l-m))/prod(1:(l+m));
    m=-m;
end
Lmn=legendre(l,cos(theta));
if l~=0
  Lmn=squeeze(Lmn(m+1,:))';
end
m=oldm;
Y=coeff*sqrt((2*l+1)*prod(1:(l-m))/(prod(1:(l+m))*4*pi)) ...
  *Lmn.*exp(i*m*phi);
