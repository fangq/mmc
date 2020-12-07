function result=dcs_g1_Db_fms(Db,tau,g1,sdsep,mu_a,mu_sp,alpha)
% 
% this function generates analytical auto-correlation decay for a semi-
% infinite medium, and computes the rms error for fminsearch fitting
%
% sdsep, mu_a and mu_sp units must be mm

zo=(mu_sp+mu_a)^-1; zb= 1.76/mu_sp; % this assumes n_tissue=1.35
r1=(sdsep^2+zo^2)^(1/2); r2=(sdsep^2+(zo+2*zb)^2)^(1/2); 
k0=2*pi*1.35/(785*1e-6); %wavelength=786e-6 mm!
kd=(3*mu_sp*mu_a+mu_sp^2*k0^2*alpha.*(6*Db.*tau)).^(1/2);
kd0=(3*mu_sp*mu_a).^(1/2);
fit_g1=exp(-kd.*r1)/r1-exp(-kd.*r2)/r2;
fit_g1_norm=exp(-kd0.*r1)/r1-exp(-kd0.*r2)/r2;

result=sum((g1-fit_g1/fit_g1_norm).^2);