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

addpath('../../matlab/');

% define 1 source and 4 detectors for a uniform cube

clear cfg newcfg

cfg.nphoton=1e6;
cfg.seed=1648335518;
[cfg.node,cfg.elem]=genT6mesh(0:2:60,0:2:60,0:2:60);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[30.1,30.2,0];
cfg.srcdir=[0 0 1];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.prop=[0 0 1 1;0.005 1.0 0.01 1.0];
cfg.debuglevel='TP';
cfg.method='elem';
cfg.isreflect=0;
cfg.detpos=[30. 20. 0. 1.
   30. 40. 0. 1.
   20. 30. 0. 1.
   40. 30. 0. 1.];
cfg.issaveexit=1;  % save detected photon exit position and angles
cfg.issaveseed=1;  % save detected photon seeds to replay later
% cfg.method=0;

newcfg=mmclab(cfg,'prep');  % preprocessing of the mesh to get the missing fields

[cube, detp, ncfg, seeds]=mmclab(newcfg); % initial simulation

% set up for wl replay

newcfg.replaydet=1;     % replay photons detected by det#1
newcfg.seed=seeds.data(:,detp.data(1,:)==newcfg.replaydet);
detp.ppath=detp.ppath(detp.data(1,:)==newcfg.replaydet,:);
detp.w0=detp.w0(detp.data(1,:)==newcfg.replaydet,:);
detp.data=detp.data(:,detp.data(1,:)==newcfg.replaydet);
% calculate the detected photon weight using the partial path output and prop
newcfg.replayweight=mmcdetweight(detp,newcfg.prop);
newcfg.replaytime=mmcdettime(detp,newcfg.prop);
newcfg.isnormalized=0;
newcfg.outputtype='wl';    % replay and get wl

% now replay the detected photons

[cube2,detp2]=mmclab(newcfg);

% the two detected photon arrays should be the same. however, because
% the program uses multi-threading, the orders may be different

if(all(ismember(round(detp.data'*1e10)*1e-10,round(detp2.data'*1e10)*1e-10,'rows')))
   disp('replay is successful :-)');
else
   disp('replay failed :-(');
end

% plot the Jacobian 3D profiles
figure;
plotmesh([cfg.node,log10(cube2.data)],cfg.elem(:,1:4),'z=0.2','linestyle','none');
hold on;
plotmesh([cfg.node,log10(cube2.data)],cfg.elem(:,1:4),'x=30','linestyle','none');
view(3)
set(gca,'xlim',[0 60])
set(gca,'ylim',[0 60])
set(gca,'zlim',[0 60])
shading interp

% set up for wp replay

newcfg.outputtype='wp';    % replay and get wp
[cube3,detp3,~,~]=mmclab(newcfg);

% the two detected photon arrays should be the same. however, because
% the program uses multi-threading, the orders may be different

if(all(ismember(round(detp.data'*1e10)*1e-10,round(detp3.data'*1e10)*1e-10,'rows')))
   disp('replay is successful :-)');
else
   disp('replay failed :-(');
end

% plot the Jacobian 3D profiles
figure;
plotmesh([cfg.node,log10(cube3.data)],cfg.elem(:,1:4),'z=0.2','linestyle','none');
hold on;
plotmesh([cfg.node,log10(cube3.data)],cfg.elem(:,1:4),'x=30','linestyle','none');
view(3)
set(gca,'xlim',[0 60])
set(gca,'ylim',[0 60])
set(gca,'zlim',[0 60])
shading interp