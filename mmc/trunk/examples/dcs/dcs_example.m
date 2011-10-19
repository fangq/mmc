addpath('../../matlab');

Db=1e-7; % simulated brownian motion diffusion coefficient
sim_beta=0.4;

mua=0.01;
mus=1;
g=0.01;
musp=mus*(1-g);
sdsep=15;
prob_scat_moving_scat=1;


[tau,g1]=generate_g1('dcs.mch',logspace(-8,0,200),'brownian',Db);

figure(1);
for I=1:4,
    subplot(2,2,I); 
    semilogx(tau,g1(I,:),'.');
    x=fminsearch(@(x) dcs_g1_Db_fms(x,tau,g1(I,:),sdsep,mua,musp,prob_scat_moving_scat),[1e-5]);
    hold on;
    semilogx(tau,dcs_g1_Db(x,tau,sdsep,mua,musp,prob_scat_moving_scat),'r');
    xlabel('\tau (s)'); ylabel('g_1(\tau)');axis tight;
    title(['Det ' num2str(I) ': fit Db: ' sprintf('%.2g',x) ' err: ' num2str((x-Db)/Db*100,'%.1f') '%']);
end

figure(2); 
for I=1:4,
    subplot(2,2,I); 
    g2=1+sim_beta*g1(I,:).^2;
    semilogx(tau,g2,'.');
    x=fminsearch(@(x) dcs_g2_Db_fms(x,tau,g2,sdsep,mua,musp,prob_scat_moving_scat),[1e-5 0.5]);
    hold on;
    semilogx(tau,dcs_g2_Db(x(1),tau,sdsep,mua,musp,prob_scat_moving_scat,x(2)),'r');
    xlabel('\tau (s)'); ylabel('g_2(\tau)');axis tight;
    title(['Det ' num2str(I) ': fit Db: ' sprintf('%.2g',x(1)) ' err: ' num2str((x(1)-Db)/Db*100,'%.1f') '%']);
end

