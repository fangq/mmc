clear;

addpath('../../matlab/');

%% testing albedo method

clear cfg newcfg

cfg.nphoton=1e7;
cfg.seed=-1;

[cfg.node,cfg.elem]=genT6mesh(0:2:60,0:2:60,0:2:60);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1,30.2,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.prop=[0 0 1 1;0.1 1.0 0.0 1.0];
cfg.debuglevel='TP';
cfg.isreflect=0;
cfg.detpos=[30.1 25.2 0. 1.0];
cfg.issaveexit=1;  % save detected photon exit position and angles
cfg.issaveseed=1;  % save detected photon seeds to replay later
cfg.isnormalized=1;
cfg.method='elem';

% albedo method
cfg.mcmethod=1; 
newcfg=mmclab(cfg,'prep');  % preprocessing of the mesh to get the missing fields
[cube1, detp1]=mmclab(newcfg); % initial simulation

% mBLL method
cfg.mcmethod=0; 
newcfg=mmclab(cfg,'prep');  % preprocessing of the mesh to get the missing fields
[cube2, detp2]=mmclab(newcfg); % initial simulation

%% CW contour comparison

figure(1);
[xi,yi] = meshgrid(0:60,0:60);
slice = 30;

node = cfg.node;
elem = cfg.elem;

[cutpos,cutvalue,facedata] = qmeshcut(elem(:,1:4),node,cube1.data,[slice 0 0; slice 0 1; slice 1 0]);
vi=griddata(cutpos(:,2),cutpos(:,3),cutvalue,xi,yi);

[cutpos,cutvalue,facedata] = qmeshcut(elem(:,1:4),node,cube2.data,[slice 0 0; slice 0 1; slice 1 0]);
vi2=griddata(cutpos(:,2),cutpos(:,3),cutvalue,xi,yi);

contour(log10(abs(vi)),[1:0.5:8],'r:','LineWidth',1.5)
hold on
contour(log10(abs(vi2)),[1:0.5:8],'b','LineWidth',1.5)
grid on

xlabel('y(mm)','FontSize',20);
ylabel('z(mm)','FontSize',20);
title('plane x=30 mm','FontSize',20);
set(gca,'FontSize',20);

legend('albedo','mBLL');

%% detected weight comparison CW

mua=0.1;
mus=1;
mut=mua+mus;

w1 = detp1.w0.*(mus/mut).^double(detp1.nscat);
w2 = detp2.w0.*exp(-mua*double(detp2.ppath));

abs(sum(w1)-sum(w2))/sum(w2)

%% detected weight comparison (TD)

c0 = 3e11;
n = 1;
cc = c0/n;

nTG = 50;
gate_width = 2e-11;

tpsf1 = zeros(nTG,1);
tpsf2 = zeros(nTG,1);

for i = 1:length(detp1.ppath)
    tg = floor(detp1.ppath(i)/cc/gate_width)+1;
    tpsf1(tg) = tpsf1(tg)+w1(i);
end

for i = 1:length(detp2.ppath)
    tg = floor(detp2.ppath(i)/cc/gate_width)+1;
    tpsf2(tg) = tpsf2(tg)+w2(i);
end

%% weighted ppath and nscat comparison

% scattering events
p1 = sum(double(detp1.nscat).*w1);
p1_ave = p1/length(detp1.nscat);
p2 = sum(double(detp2.nscat).*w2);
p2_ave = p2/length(detp2.nscat);

% pathlength
L1 = sum(double(detp1.ppath).*w1);
L1_ave = L1/length(detp1.nscat);
L2 = sum(double(detp2.ppath).*w2);
L2_ave = L2/length(detp2.nscat);
