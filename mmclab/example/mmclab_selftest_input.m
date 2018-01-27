%%-----------------------------------------------------------------
%% MMClab self-test script
%%-----------------------------------------------------------------
%
% In this example, we generate various true or false input conditions
% and test the robustness of the script to these conditions 
%

dat=mmclab();

cfg=struct();
dat=mmclab(cfg);

[cfg.node face cfg.elem]=meshabox([0 0 0], [10 10 10],1);
dat=mmclab(cfg);

cfg.elem(:,5)=1;
dat=mmclab(cfg);

cfg.srcpos=[5 5 0];
cfg.srcdir=[0 0 1];
cfg.nphoton=1000;
cfg.prop=[0 0 1 1;0.1 1 0.2 1.3];
dat=mmclab(cfg);

cfg.tstart=0;
cfg.tend=1e-9;
dat=mmclab(cfg);

cfg.tstep=10e-9;
dat=mmclab(cfg);

cfg.srcpos=[-1 -1 -1];
dat=mmclab(cfg);

cfg.srcparam1=[0 0 0 0];
cfg.srcparam2=[0 0 0 0];
dat=mmclab(cfg);
