%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input

cfg.nphoton=1e6;
[cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],10);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1 30.1 0];
cfg.srcdir=[0 0 1];
cfg.prop=[0 0 1 1;0.005 1 0 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';
cfg.basisorder=2;
cfg.isatomic=1; % set to 0 to disable atomic to gain speed
%cfg.method='P';


% run the simulation

tic;
[flux, detps, newcfg]=mmclab(cfg);
toc;

% plotting the result

plotmesh([cfg.node(:,1:3),log10(abs(flux.data(1:size(cfg.node,1))))],cfg.elem,'y=30','facecolor','interp','linestyle','none')
view([0 1 0]);
colorbar;
