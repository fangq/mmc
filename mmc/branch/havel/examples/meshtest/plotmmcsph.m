%%-----------------------------------------------------------------
%% add paths to the necessary toolboxes
%%-----------------------------------------------------------------

addpath('/space/kwafoo/2/users/fangq/Projects/mcx/utils/')
addpath('../../matlab')

c0=299792458000;
twin=[5e-11:1e-10:5e-9];
gates=50;
clines=[-1:0.5:8];

[xi,yi]=meshgrid(1:60,0:60);

%%-----------------------------------------------------------------
%% load MCX results
%%-----------------------------------------------------------------
mcx=loadmc2('../mcxsph/spherebox.mc2', [60 60 60 gates]);
cwmcx=sum(mcx,4);

%%-----------------------------------------------------------------
%% generate/load analytical solution for sphere inside infinite slab
%%-----------------------------------------------------------------

%[phi_ana,xa,ya,za]=sphdiffusionslab(0,0,60,-22:0.8:22,0,-30:0.8:10);
%save sphdiffsemiinf.mat phi_ana xa ya za

load sphdiffsemiinf.mat
idx=find((xa(:)<-12 | xa(:)>12) & za(:)>-5);
phi_ana(idx)=nan;
idx=find((xa(:)<-10 | xa(:)>10) & za(:)>0);
phi_ana(idx)=nan;

%%-----------------------------------------------------------------
%% generate the contour of the inclusion
%%-----------------------------------------------------------------

[xcirc,ycirc] = cylinder(10,200);
xcirc=xcirc(1,:)+30;
ycirc=ycirc(1,:)+30;

% create the voxel-contour of the sphere for MCX
%dim=60;
%[xv,yv,zv]=meshgrid(1:dim,1:dim,1:dim);
%dist=(xv-30).^2+(yv-30).^2+(zv-30).^2;

%vv=zeros(size(xv));
%vv(dist<100)=1;

%c=contourc(squeeze(vv(:,30,:)),1);
%plot(c(1,2:end),c(2,2:end),'c--')

%%-----------------------------------------------------------------
%% plot sphere 1
%%-----------------------------------------------------------------
node1=readmmcnode('node_sph1.dat');  
elem1=readmmcelem('elem_sph1.dat');
load sph1.dat
sph1=reshape(sph1(:,end),[size(node1,1),length(sph1)/size(node1,1)]);
s1=sum(sph1,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem1(:,1:4),node1,s1,[0 30 0; 0 30 1; 1 30 0]);
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])

vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana))+10,clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(squeeze(abs(cwmcx(:,30,:)))'),clines,'b-')
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MCX','MMCM')
legend boxoff;
box on;
set(gcf,'PaperPositionMode','auto');
print -depsc2 sph1.eps

%%-----------------------------------------------------------------
%% plot sphere 2
%%-----------------------------------------------------------------

node2=readmmcnode('node_sph2.dat');						    
elem2=readmmcelem('elem_sph2.dat');										    
load sph2.dat
sph2=reshape(sph2(:,end),[size(node2,1),length(sph2)/size(node2,1)]);
s2=sum(sph2,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem2(:,1:4),node2,s2,[0 30 0; 0 30 1; 1 30 0]);
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])
vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana))+10,clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(squeeze(abs(cwmcx(:,30,:)))'),clines,'b-')
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MCX','MMCM')
legend boxoff;
box on;
set(gcf,'PaperPositionMode','auto');
print -depsc2 sph2.eps

%%-----------------------------------------------------------------
%% plot sphere 3
%%-----------------------------------------------------------------

node3=readmmcnode('node_sph3.dat');
elem3=readmmcelem('elem_sph3.dat');
load sph3.dat
sph3=reshape(sph3(:,end),[size(node3,1),length(sph3)/size(node3,1)]);
s3=sum(sph3,2);
[cutpos,cutvalue,facedata]=qmeshcut(elem3(:,1:4),node3,s3,[0 29 0; 0 29 1; 1 29 0]);

vi=griddata(cutpos(:,1),cutpos(:,3),cutvalue,xi,yi);

figure
%patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',log10(cutvalue),'FaceColor','interp','linestyle','none');
%view([0 1 0])
hold on
[cc,hc]=contour(xa+30,za+31,log10(abs(phi_ana))+10,clines,'color',[0.7 0.7 0.7],'linewidth',2);
contour(log10(squeeze(abs(cwmcx(:,30,:)))'),clines,'b-')
contour(log10(abs(vi)),clines,'r:')
plot(xcirc,ycirc,'k--','linewidth',2);

axis equal
set(gca,'xlim',[1 60]);
set(gca,'ylim',[1 60]);
set(gca,'fontsize',18)
xlabel('x (mm)')
ylabel('z (mm)')
legend('Diffusion','MCX','MMCM')
legend boxoff;
box on;

set(gcf,'PaperPositionMode','auto');
print -depsc2 sph3.eps

