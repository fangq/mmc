%%-----------------------------------------------------------------
% In this example, we demonstrate light transport simulation in a mesh-based 
% head atlas generated from the USC 19.5-year group atlas template[Sanchez2012]). 
%
% This script is a lightweight version of the DMMC[Yan2019] simulation demo used 
% for Fig.9(a) in TranYan2019(submitted). For simplicity, a coarse mesh (lower  
% sampling rate) that represents part of the head(z>155) is simulated.
% 
% [Sanchez2012] C.E.Sanchez J.E.Richards and C.R.Almli, “Age-Specific MRI Templates
% for Pediatric Neuroimaging,” Developmental Neuropsychology 37, 379–399 (2012).
%
% [Yan2019] S.Yan, A.P.Tran, and Q.Fang, “Dual-grid mesh-based monte carlo algorithm  
% for efficient photon transport simulations in complex three-dimensional media,” 
% Journal of Biomedical Optics 24(2), 020503 (2019).
%
% [TranYan2019] A.P.Tran, S.Yan and Q.Fang, "Improving model-based fNIRS 
% analysis using mesh-based anatomical and light-transport models".
%
%%-----------------------------------------------------------------

clc
clear
load('head_atlas.mat');
%% prepare cfg for MMCLAB simulation
clear cfg
cfg.nphoton=1e7; %takes ~240s on Intel i7-8700K(12 threads)

% medium labels:0-ambient air,1-air cavities,2-scalp,3-skull,4-csf,5-gray matter,6-white matter 
cfg.node=double(node);
cfg.elem=double(elem);
cfg.elemprop=double(prop);
cfg.prop=[0,0,1,1;0,0,1,1;0.019 7.8 0.89 1.37;0.019 7.8 0.89 1.37;0.004 0.009 0.89 1.37;0.02 9.0 0.89 1.37;0.08 40.9 0.84 1.37];

% light source
cfg.srctype='pencil';
cfg.srcdir=[-0.5086,-0.1822,-0.8415]; %inward-pointing source
cfg.srcpos=[133.5370,90.1988,200.0700]; %pencil beam source placed at EEG 10-5 landmark:"C4h"
    cfg.srcpos=cfg.srcpos+0.001*cfg.srcdir; %ensure source is inside the mesh

% time windows
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;

% other simulation parameters
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.method='grid'; %DMMC mode

%% run mmc simulation
[flux]=mmclab(cfg);

%% post-simulation data processing and visualization
% convert time-resolved fluence to CW fluence
CWfluence=sum(flux.data*cfg.tstep,4); 

% coronal plane selected for fluence plot
y_plane=90.5;   

% fix spatial offset between the coordinate origins of DMMC AABB and the mesh
DMMC_AABB_origin=min(cfg.node); %spatial offset of the DMMC AABB
DMMC_voxel_origin=DMMC_AABB_origin+0.5; %centroid of the 1st voxel for DMMC AABB
DMMC_AABB_dimension=size(CWfluence);
CWfluence_DMMC=squeeze(CWfluence(:,ceil(y_plane-DMMC_AABB_origin(1,2)),:)); %along y-axis, DMMC AABB happens to have an offset very close to 3(2.9955)
[xx,zz]=meshgrid(DMMC_voxel_origin(1,1):(DMMC_voxel_origin(1,1)+DMMC_AABB_dimension(1,1)-1),DMMC_voxel_origin(1,3):(DMMC_voxel_origin(1,3)+DMMC_AABB_dimension(1,3)-1));

figure;
clines=-20:0.5:0;
contourf(xx,zz,log10(abs(CWfluence_DMMC')),clines,'linestyle','--','color','w','linewidth',1.5,'DisplayName','DMMC');
hold on;axis equal;
colorbar('EastOutside');

% plot tissue boundaries, auxiliary functions are part of iso2mesh(http://iso2mesh.sf.net)
plane=[0 y_plane 0; 0 y_plane 1; 1 y_plane 0];
label=unique(prop);
face=[];
for i=1:2:size(label,1)
    newface=volface(cfg.elem(prop==label(i),:));
    face=[face;newface];
end
[cutpos,cutvalue,cutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
[cutpos,cutedges]=removedupnodes(cutpos,cutedges);
cutloop=extractloops(cutedges);
[nanidx]=find(isnan(cutloop));

for i=1:size(nanidx,2)
    if(i==9) 
        continue;
    end
    if(i==1)
        plot(cutpos(cutloop(1:(nanidx(i)-1)),1),cutpos(cutloop(1:(nanidx(i)-1)),3),'linestyle','--','color','k','LineWidth',1.25,'HandleVisibility','off');
    elseif(i==2)
        plot(cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),1),cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),3),'linestyle','--','color','k','LineWidth',1.25,'HandleVisibility','off');
    else
        plot(cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),1),cutpos(cutloop((nanidx(i-1)+1):(nanidx(i)-1)),3),'linestyle','--','color','k','LineWidth',1.25,'HandleVisibility','off');
    end
end

% plot source, legend, etc.
plot(cfg.srcpos(1,1),cfg.srcpos(1,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'DisplayName','source');
lg=legend('Location','northeast');
set(lg,'color','[0.85 0.85 0.85]');
set(lg,'box','on');
set(gca,'ylim', [160 225]);ylabel('z [mm]');
set(gca,'xlim', [10 165]);xlabel('x [mm]');
set(gca,'clim',[-12 0]);
set(gca,'fontsize',18);
set(gca, 'FontName', 'Times New Roman');