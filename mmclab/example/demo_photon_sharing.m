% In this example, we show the photon sharing feature of MMCLAB.
% i.e. obtain results of multiple source patterns with only one simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare simulation input

cfg.nphoton=3e6;
[cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 20],2, 2);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[10 10 -2];
cfg.srcdir=[0 0 1];
cfg.srctype='pattern';
cfg.srcparam1=[40.0 0.0 0.0 40];
cfg.srcparam2=[0.0 40.0 0.0 40];
cfg.prop=[0 0 1 1;0.01 10 0.9 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.debuglevel='TP';

%% 5 patterns represented by a 3D matrix

pat=zeros(40,40,5);
pat(:,:,1) = ones(40);  % full field illumination pattern
pat(1:20,:,2) = 1;      % pattern with left half bright
pat(:,1:20,3) = 1;      % pattern with top half bright
pat(11:30,11:30,4) = 1; % pattern with bright square in the middle
pat(16:25,:,5) = 1;  
pat(:,16:25,5) = 1;     % pattern with a bright cross

cfg.srcpattern = pat;

%% run the simulation

flux=mmclab(cfg);

%% plot results (same as in the mmc example)

[node, ~, elem]=meshabox([0 0 0],[60 60 20],2, 2);
data = flux.data(:,1:length(node),:);
cwdata = squeeze(sum(data,3));

figure();

subplot(1,3,1); title('pattern 1');
plotmesh([node,(cwdata(1,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,2); title('pattern 2');
plotmesh([node,(cwdata(2,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,3); title('pattern 3');
plotmesh([node,(cwdata(3,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,5); title('pattern 4');
plotmesh([node,(cwdata(4,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,6); title('pattern 5');
plotmesh([node,(cwdata(5,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])
