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
    0.0190    7.8182    0.8900    1.3700 % scalp
    0.0190    7.8182    0.8900    1.3700 % skull
    0.0040    0.0090    0.8900    1.3700 % csf
    0.0200    9.0000    0.8900    1.3700 % gray matters
    0.0800   40.9000    0.8400    1.3700 % white matters
         0         0    1.0000    1.0000]; % air pockets

cfg.seed=29012392;
cfg.nphoton=1e7;
cfg.method='grid';

%%
cfg.gpuid=1;

b3gpu=mmclab(cfg);
b3gpu=sum(b3gpu.data,4)*cfg.tstep;
%%
cfg.gpuid=-1;

b3cpu=mmclab(cfg);
b3cpu=sum(b3cpu.data,4)*cfg.tstep;


%%-----------------------------------------------------------------
%% generate a contour plot along y=30.2
%%-----------------------------------------------------------------
figure

clines = 8:-0.5:-8;

hold on
[c h2]=contourf(log10(squeeze(b3cpu(75,:,:))'),clines,'k-');
contour(log10(abs(squeeze(b3gpu(75,:,:)))'),clines,'y--')

axis equal  
set(gca,'xlim',[1 182])
set(gca,'fontsize',20)
xlabel('x (mm)')
ylabel('z (mm)')
legend('MMC','MMCL')
legend boxoff;
box on;

