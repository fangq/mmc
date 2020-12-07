clear all;

%% 2-layer mesh model

layercount=2;

if(layercount==2)
    [node,face,c0]=latticegrid([0 60],[0 60],[0 25 30]);
    c0(:,4)=[2;1];   % maximum element size for bottom (label 1) and top (label 2) layers
    
    % simulate the optical properties of skull and gray-matter of a human brain model
    % see http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
    % cfg.prop=[0 0 1 1;0.02 9.0, 0.89 1.37;0.019 7.8 0.89 1.37]; 
else
    [node,face,c0]=latticegrid([0 60],[0 60],[0 20 25 30]); % if you like a 3-layer model
    c0(:,4)=[2;2;1];
    % cfg.prop=[0 0 1 1;0.02 9.0, 0.89 1.37;0.004 0.009, 0.89 1.37;0.019 7.8 0.89 1.37]; 
end
[no,el]=surf2mesh(node,face,[],[],1,[],c0);

figure(1); plotmesh(no,el);
savemmcmesh('media',no,el);

%% add source and retessellate mesh

cfg.srctype='fourier';
cfg.srcpos=[10 10 35];
kx=3;                       % wave number in the x-dir
ky=0;                       % wave number in the x-dir
xphase=pi/3;                % phase offset in the x-dir, must < 2pi
yphase=0;                   % phase offset in the x-dir, must < 2pi
cfg.srcparam1=[40 0 0 kx+xphase/(2*pi)];   % 3 is k-number in x direction
cfg.srcparam2=[0 40 0 ky+yphase/(2*pi)];
cfg.srcdir=[0 0 -1];

srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
    'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);

[node,elem] = mmcaddsrc(no,el,mmcsrcdomain(srcdef,[min(no);max(no)]));

figure(2); plotmesh(node,elem);

% save the final re-tessellated mesh
savemmcmesh('sfdi',node,elem);
