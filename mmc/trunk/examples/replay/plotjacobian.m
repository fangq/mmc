% verify the replay outcome and plot the Jacobian profile

% the photons detected in the original run
[data1,header1]=loadmch('step1.mch');
header1

% the photons detected in the wl replay phase
[data2,header2]=loadmch('step2.mch');
header2

% the photons detected in the wp replay phase
[data3,header3]=loadmch('step3.mch');
header3

data1=data1(data1(:,1)==1,:);

% the two sets of captured photons should be identical
if(all(ismember(round(data1*1e10)*1e-10,round(data2*1e10)*1e-10,'rows')) && ...
   all(ismember(round(data1*1e10)*1e-10,round(data2*1e10)*1e-10,'rows')))
   disp('replay is successful :-)');
else
   disp('replay failed :-(');
end

% plot the wl profile

load step2.dat
[no,el]=readmmcmesh('replaytest');

figure;
plotmesh([no,log10(step2(:,2))],el(:,1:4),'z=0.2','linestyle','none');
hold on;
plotmesh([no,log10(step2(:,2))],el(:,1:4),'x=30','linestyle','none');
view(3)
set(gca,'xlim',[0 60])
set(gca,'ylim',[0 60])
set(gca,'zlim',[0 60])
shading interp

% plot the wp profile

load step3.dat
[no,el]=readmmcmesh('replaytest');

figure;
plotmesh([no,log10(step3(:,2))],el(:,1:4),'z=0.2','linestyle','none');
hold on;
plotmesh([no,log10(step3(:,2))],el(:,1:4),'x=30','linestyle','none');
view(3)
set(gca,'xlim',[0 60])
set(gca,'ylim',[0 60])
set(gca,'zlim',[0 60])
shading interp