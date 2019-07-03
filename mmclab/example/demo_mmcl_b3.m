addpath('../../../matlab/'); 

if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.tar.gz','file'))
          urlwrite('http://downloads.sourceforge.net/project/mcx/mmc/AdultBrain%20Atlas%20FEM%20Mesh/Version%202/MMC_Collins_Atlas_Mesh_Version_2L.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fmcx%2Ffiles%2Fmmc%2FAdultBrain%2520Atlas%2520FEM%2520Mesh%2FVersion%25202%2F&ts=1451344236&use_mirror=tcpdiag','MMC_Collins_Atlas_Mesh_Version_2L.tar.gz');
     end
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.tar','file'))
          gunzip('MMC_Collins_Atlas_Mesh_Version_2L.tar.gz');
     end
     if(~exist('MMC_Collins_Atlas_Mesh_Version_2L.mat','file'))
          untar('MMC_Collins_Atlas_Mesh_Version_2L.tar');
     end
end

load MMC_Collins_Atlas_Mesh_Version_2L.mat

clear cfg
cfg.node=node;
cfg.elem=elem(:,1:4);
cfg.elemprop=elem(:,5);

cfg.tstart=0;
cfg.tend=5e-09;
cfg.tstep=2e-10;

cfg.srcpos=[75 67.38 167.5];
cfg.srcdir=[0.1636 0.4569 -0.8743];
cfg.srcdir=cfg.srcdir/norm(cfg.srcdir);

cfg.detpos=[   75.0000   77.1900  170.3000    3.0000
   75.0000   89.0000  171.6000    3.0000
   75.0000   97.6700  172.4000    3.0000
   75.0000  102.4000  172.0000    3.0000];

cfg.prop=[         0         0    1.0000    1.0000 % background/air
    0.0190    7.8182    0.8900    1.3700 % scalp & skull are both labeled using 1
%     0.0190    7.8182    0.8900    1.3700 % skull
    0.0040    0.0090    0.8900    1.3700 % csf
    0.0200    9.0000    0.8900    1.3700 % gray matters
    0.0800   40.9000    0.8400    1.3700 % white matters
         0         0    1.0000    1.0000]; % air pockets

cfg.seed=29012392;
cfg.nphoton=1e8;
cfg.method='elem';

%%
cfg.gpuid=1;

b3gpu=mmclab(cfg);
b3gpu=sum(b3gpu.data,2)*cfg.tstep;
%%
cfg.gpuid=-1;

b3cpu=mmclab(cfg);
b3cpu=sum(b3cpu.data,2)*cfg.tstep;

%%-----------------------------------------------------------------
%% generate a contour plot along y=74.5
%%-----------------------------------------------------------------
figure

clines = 0:-0.5:-10;

x_min=min(cfg.node(:,1)); x_max=max(cfg.node(:,1)); bb_size(1,1)=floor(x_max-x_min);
y_min=min(cfg.node(:,2)); y_max=max(cfg.node(:,2)); bb_size(1,2)=floor(y_max-y_min);
z_min=min(cfg.node(:,3)); z_max=max(cfg.node(:,3)); bb_size(1,3)=floor(z_max-z_min);
[yi,zi]=ndgrid((y_min+0.5):(y_min+0.5+bb_size(1,2)-1),(z_min+0.5):(z_min+0.5+bb_size(1,3)-1));

[cutpos1,cutvalue1,facedata1]=qmeshcut(elem(:,1:4),node,b3cpu,[74.5 0 0;74.5 0 1;74.5 1 0]);
[cutpos2,cutvalue2,facedata2]=qmeshcut(elem(:,1:4),node,b3gpu,[74.5 0 0;74.5 0 1;74.5 1 0]);
b3cpu_interpolated=griddata(cutpos1(:,2),cutpos1(:,3),cutvalue1,yi,zi);
b3gpu_interpolated=griddata(cutpos2(:,2),cutpos2(:,3),cutvalue2,yi,zi);

b3cpu_interpolated(find(b3cpu_interpolated==0))=nan;
b3gpu_interpolated(find(b3gpu_interpolated==0))=nan;

hold on
contourf(yi,zi,log10(b3cpu_interpolated),clines,'k-','linewidth',2);
contour(yi,zi,log10(b3gpu_interpolated),clines,'g--','linewidth',2)

x_cut=75-0.5;
plane=[x_cut 0 0;x_cut 0 1;x_cut 1 0];
[cutpos,cutvalue,cutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
[cutpos,cutedges]=removedupnodes(cutpos,cutedges);
cutloop=extractloops(cutedges);
[nanidx]=find(isnan(cutloop));

for i=1:size(nanidx,2)
    if(i==1)
        plot(cutpos(cutloop(1:(nanidx(i)-1)),2),cutpos(cutloop(1:(nanidx(i)-1)),3),'color','k','LineWidth',1.25,'HandleVisibility','off');
    elseif(i==2)
        plot(cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),2),cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),3),'linestyle','--','color','k','LineWidth',1.25,'HandleVisibility','off');
    else
        plot(cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),2),cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),3),'linestyle','--','color','k','LineWidth',1.25,'HandleVisibility','off');
    end
end

axis equal  
set(gca,'xlim',[3 217]);
set(gca,'ylim',[-1 175]);
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('MMC','MMCL')
legend boxoff;
box on;
