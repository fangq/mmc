addpath('../../matlab/');
addpath('..');

clear cfg

[node,face,c0]=latticegrid([0 60],[0 60],[0 5 10]);
c0(:,4)=[2;3];   % maximum element size for bottom (label 1) and top (label 2) layers
[node,elem]=surf2mesh(node,face,[],[],1,[],c0);

detdef=struct('srctype','planar','srcpos',[10,10,-1],'srcdir',[0 0 1],...
                'srcparam1',[40 0 0 40],'srcparam2',[0 40 0 40]);
            
[cfg.node,cfg.elem]=mmcadddet(node,elem,detdef);
cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);
cfg.prop = [0 0 1 1;
                0.01 1.0 0.01 1.37
                0.05 10.0 0.9 1.37];
    
cfg.srcpos = [25.0 35.0 10.0];
cfg.e0 = tsearchn(cfg.node,cfg.elem(:,1:4),cfg.srcpos);
cfg.srcdir = [0, 0, -1];

cfg.srctype = 'pencil';
cfg.srcparam1 = [0 0 0 0];
cfg.srcparam2 = [0 0 0 0];

cfg.detpos = [10.0,10.0,-1.0,0];
cfg.detparam1 = [40.0 0.0 0.0 40];
cfg.detparam2 = [0.0 40.0 0.0 40];

cfg.tstart = 0;
cfg.tend = 2e-9;
cfg.tstep = 4e-11;

cfg.nphoton = 1e6;
cfg.seed = 12345678;
cfg.debuglevel = 'TP';
cfg.issaveexit = 2;

[flux,detp]=mmclab(cfg);

figure;
imagesc(sum(detp.data,3)');
axis equal;
