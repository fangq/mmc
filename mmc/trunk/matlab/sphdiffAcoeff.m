function A=sphdiffAcoeff(m,l,cfg)
if((cfg.src(2)==pi | cfg.src(2)==0) & m~=0)
    A=0;
    return;
end
x=cfg.kout*cfg.a;
y=cfg.kin*cfg.a;
Dout=cfg.Dout;
Din=cfg.Din;
A=-i*cfg.v*cfg.kout/Dout*spbesselh(l,1,cfg.kout*cfg.src(1))*conj(spharmonic(l,m,cfg.src(2),cfg.src(3)))...
  *(Dout*x*spbesseljprime(l,x)*spbesselj(l,y)-Din*y*spbesselj(l,x)*spbesseljprime(l,y))...
  /(Dout*x*spbesselhprime(l,1,x)*spbesselj(l,y)-Din*y*spbesselh(l,1,x)*spbesseljprime(l,y));
