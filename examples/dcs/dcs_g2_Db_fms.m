function result=dcs_g2_Db_fms(s,tau,g2,sdsep,mu_a,mu_sp,alpha)
% 
% this function generates analytical field auto-correlation decay for a semi-
% infinite medium, then uses the Siegert relationship to derive the
% intensity temporal auto-correlation decay assuming ergodicity
%
% sdsep, mu_a and mu_sp units must be mm

Db=s(1);
beta=s(2);

zo=(mu_sp+mu_a)^-1; zb= 1.76/mu_sp; % this assumes n_tissue=1.35
r1=(sdsep^2+zo^2)^(1/2); r2=(sdsep^2+(zo+2*zb)^2)^(1/2); 
k0=2*pi*1.35/(785*1e-6); %wavelength=786e-6 mm!
kd=(3*mu_sp*mu_a+mu_sp^2*k0^2*alpha.*(6*Db.*tau)).^(1/2);
kd0=(3*mu_sp*mu_a).^(1/2);
fit_g1=exp(-kd.*r1)/r1-exp(-kd.*r2)/r2;
fit_g1_norm=exp(-kd0.*r1)/r1-exp(-kd0.*r2)/r2;

g1_norm=fit_g1/fit_g1_norm;
fit_g2=1+beta*(g1_norm).^2;

result=sum((g2-fit_g2).^2);
