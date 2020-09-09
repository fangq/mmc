%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input

cfg.nphoton=1e7;
[cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[20 30 -10];
cfg.srcdir=[0 0.2 sqrt(1-0.04)];
cfg.prop=[0 0 1 1;0.005 0.1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.unitinmm=1;
cfg.method='elem';
cfg.srctype='slit';
cfg.srcparam1=[30 0 0 0];
cfg.srcparam2=[0 0 0 0];

% run the simulation

flux=mmclab(cfg);

% plotting the result

% if you have the SVN version of iso2mesh, use the next line to plot:
% qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30','linestyle','none');

plotmesh([cfg.node(:,1:3),log10(abs(flux.data(1:size(cfg.node,1))))],cfg.elem,'x=40','facecolor','interp','linestyle','none')
hold on;
plotmesh([cfg.node(:,1:3),log10(abs(flux.data(1:size(cfg.node,1))))],cfg.elem,'y=30','facecolor','interp','linestyle','none')
view(3);
colorbar;
