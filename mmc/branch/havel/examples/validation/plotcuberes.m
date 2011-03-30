%%-----------------------------------------------------------------
%% add paths to the necessary toolboxes
%%-----------------------------------------------------------------

%addpath('/Please/add/path/to/mcx/utils/')
addpath('../../matlab')

c0=299792458000;
t0=0;
dt=1e-10;
t1=5e-9;

twin=[t0+dt/2:dt:t1];
gates=length(twin);

%%-----------------------------------------------------------------
%% load MCX results
%%-----------------------------------------------------------------
mcx=loadmc2('../mcxsph/box.mc2', [60 60 60 gates]);
cwmcx=sum(mcx,4); 

%%-----------------------------------------------------------------
%% load MMC results
%%-----------------------------------------------------------------
node=readmmcnode('node_cube.dat');
elem=readmmcelem('elem_cube.dat');

load cube.dat
cube=reshape(cube(:,end),[size(node,1),length(cube)/size(node,1)]);
cwcb=sum(cube,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node,cwcb,[0 30.2 0; 0 30.2 1; 1 30.2 0]);

[xi,yi]=meshgrid(0.5:60,0:60);
vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

%%-----------------------------------------------------------------
%% plot the time-domain TPSF at [30 14 10]
%%-----------------------------------------------------------------

figure

srcpos=[30 30 0];
detpos=[30 14 10];

hold on
semilogy((1:gates)/10,tddiffusion(0.005, 1, c0, 0, srcpos, detpos,twin),'r');
semilogy((1:gates)/10,squeeze(mcx(detpos(1),detpos(2),detpos(3),:)),'o');
semilogy((1:gates)/10,squeeze(cube(find(node(:,1)==detpos(1) & node(:,2)==detpos(2) & node(:,3)==detpos(3)),:)),'+');

set(gca,'fontsize',20)
xlabel('t (ns)')
ylabel('Fluence TPSF (1/mm^2)')
set(gca,'yscale','log')
legend('Diffusion','MCX','MMCM')
legend boxoff;
box on;

set(gcf,'PaperPositionMode','auto');
saveas(gcf,'box_td.fig');


%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure

hold on
contour(log10(squeeze(abs(cwmcx(:,30,:)))'),[-1:0.5:8],'b-')
contour(log10(abs(vi)),[-1:0.5:8],'r:')

axis equal  
set(gca,'xlim',[1 60])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('MCX','MMCM')
legend boxoff;
box on;

set(gcf,'PaperPositionMode','auto');
saveas(gcf,'box.fig');

