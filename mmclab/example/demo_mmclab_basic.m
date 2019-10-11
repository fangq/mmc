%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input

clear cfg
cfg.nphoton=1e6;
[cfg.node, face, cfg.elem]=meshabox([0 0 0],[60 60 30],6);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30 30 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.issaveref=1;  % in addition to volumetric fluence, also save surface diffuse reflectance

%% run the simulation

flux=mmclab(cfg);

%% plotting the result

% plot the fluence in log10 scale
subplot(121);
flux_dmmc=flux.data(1:end-1,1:end-1,1:end-1);
imagesc(log10(abs(squeeze(flux_dmmc(:,31,:)))));
colorbar;

% plot the surface diffuse reflectance
if(isfield(cfg,'issaveref') && cfg.issaveref==1)
    subplot(122);
    faces=faceneighbors(cfg.elem,'rowmajor');
    hs=plotmesh(cfg.node,faces,'cdata',log10(flux.dref(:,1)),'linestyle','none');
    colorbar;
end
