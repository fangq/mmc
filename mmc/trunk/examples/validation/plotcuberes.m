%%-----------------------------------------------------------------
%% add paths to the necessary toolboxes
%%-----------------------------------------------------------------

addpath('/space/kwafoo/2/users/fangq/Projects/mcx/utils/')
addpath('../../matlab')

c0=299792458000;
twin=[5e-11:1e-10:5e-9];
gates=50;

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

srcpos=[31 31 1];
detpos=[31 15 11];

hold on
semilogy((1:50)/10,tddiffusion(0.005, 1, c0, 0, srcpos, detpos,twin),'r');
semilogy((1:50)/10,squeeze(mcx(30,14,10,:)),'o');
semilogy((1:50)/10,squeeze(cube(5038,:)),'+');

set(gca,'fontsize',18)
xlabel('t (ns)')
ylabel('Fluence TPSF (1/mm^2)')
set(gca,'yscale','log')
legend('Diffusion','MCX','MMC')
box on;

saveas(gcf,'box_td.fig');


%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure

hold on
contour(log10(squeeze(abs(cwmcx(:,31,:)))'),[-1:0.5:8],'b--')
contour(log10(abs(vi)),[-1:0.5:8],'r:')

axis equal  
set(gca,'xlim',[1 60])
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('y (mm)')
legend('MCX','MMC')
box on;

saveas(gcf,'box.fig');

