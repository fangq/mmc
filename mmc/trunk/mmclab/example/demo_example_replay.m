%%-----------------------------------------------------------------
%% demonstration of the replay feature
%%-----------------------------------------------------------------
%
% In this example, we run MMC using a homogeneous cubic domain,
% same as in demo_example_validation.m. We save the seeds of the
% detected photons, and then rerun these photons again (i.e. 
% the "replay" part). This feature allows one to identify features
% that are only related to a source/detector pair.
%
%%-----------------------------------------------------------------

%addpath('/Please/add/path/to/mcx/utils/')
addpath('../../matlab/');

cfg.nphoton=1e6;
cfg.seed=27182818;
[cfg.node,cfg.elem]=genT5mesh(0:2:60,0:2:60,0:2:60);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1,30.2,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.prop=[0 0 1 1;0.005 1.0 0.01 1.0];
cfg.debuglevel='TP';
cfg.isreflect=0;
cfg.detpos=[30. 20. 0. 1.
   30. 40. 0. 1.
   20. 30. 0. 1.
   40. 30. 0. 1.];
cfg.issaveexit=1;

newcfg=mmclab(cfg,'prep');

[cube detp ncfg seeds]=mmclab(newcfg');

% now replay the detected photons

newcfg.seed=seeds.data;
newcfg.isnormalized=0;
newcfg.outputtype=4;

[cube2 detp2 ncfg2 seeds2]=mmclab(newcfg);

% the two detected photon arrays should be the same. however, because
% the program uses multi-threading, the orders may be different

[isreplayed, mapidx]=ismember(detp.data',detp2.data','rows');

if(all(isreplayed))
   disp('replay is successful :-)');
else
   disp('replay failed :-(');
end

qmeshcut(cfg.elem(:,1:4),cfg.node,sum(cube2.data,2),'y=30.2');
view([0 1 0]);

