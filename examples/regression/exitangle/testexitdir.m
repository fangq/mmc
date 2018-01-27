dat=loadmch('tank_planar.mch');
dd=dat(:,7:9);
dl=sqrt(sum(dd.*dd,2));
% plot(dl)

find(sum(dd(:,1:2).*dd(:,1:2),2)<1e-5)

