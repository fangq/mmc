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

cfg.nphoton=1e7;
cfg.seed=27182818;
[cfg.node,cfg.elem]=genT6mesh(0:2:60,0:2:60,0:2:60);
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
cfg.issaveexit=1;  % save detected photon exit position and angles
cfg.issaveseed=1;  % save detected photon seeds to replay later

newcfg=mmclab(cfg,'prep');  % preprocessing of the mesh to get the missing fields

[cube, detp, ncfg, seeds]=mmclab(newcfg); % initial simulation

% plot the detected photon positions for verification

figure;
plotmesh(detp.data(4:6,:)');

% perturb a specific element, and get the Jacobian using the perturbation
% method

id=find(cfg.node(:,1)<32 & cfg.node(:,1)>28 & cfg.node(:,2)<32 & cfg.node(:,2)>28 & cfg.node(:,3)<12 & cfg.node(:,3)>8);
eid=find(cfg.elem==id);
newcfg.elemprop(eid(1))=2;
deltamua=1e-4;
newcfg.prop(3,:)=cfg.prop(2,:)+[deltamua 0 0 0]; % perturb the mua by a small amount

cubemua=mmclab(newcfg);  % run a new forward with the perturbed optical properties
Jtest0=(cube.data(id)-cubemua.data(id))/deltamua  % method 1 for estimating the Jacobian

newcfg.elemprop(eid(1))=1;
newcfg.prop=newcfg.prop(1:2,:);

% set up for replay

newcfg.replaydet=1;     % replay photons detected by det#1
newcfg.seed=seeds.data(:,detp.data(1,:)==newcfg.replaydet);
detp.data=detp.data(:,detp.data(1,:)==newcfg.replaydet);
% calculate the detected photon weight using the partial path output and prop
newcfg.replayweight=mmcdetweight(detp.data,newcfg.prop,detp.data(end,:));
newcfg.isnormalized=0;
newcfg.outputtype='jacobian';    % replay and get the Jacobian
det1val=mean(newcfg.replayweight);

% now replay the detected photons

[cube2 detp2 ncfg2 seeds2]=mmclab(newcfg);

Jtest1=cube2.data(id)   % method 2 for estimating the sensitivity

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

% Method 3 - use perturbation MC in the replay mode
newcfg.elemprop(eid(1))=2;
newcfg.prop(3,:)=cfg.prop(2,:)+[0.0001 0 0 0]; % perturb the mua by a small amount
newcfg.outputtype='flux';    % replay and get the Jacobian

[cube3 detp3 ncfg3 seeds3]=mmclab(newcfg);

det1val2=mean(mmcdetweight(detp3.data,newcfg.prop,detp3.data(end,:)));
Jtest2=(det1val-det1val2)/deltamua   % method 3 for estimating the sensitivity, missing a volume?


