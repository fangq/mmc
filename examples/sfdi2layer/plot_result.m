% CW fluence
data = load('sfdi.dat');
data = reshape(data(:,2),[],50);
layercw = sum(data(1:end-3,:),2)*1e-10;

% media mesh
cfg.node = readmmcnode('node_media.dat');
cfg.elem = readmmcelem('elem_media.dat');
cfg.elemprop=cfg.elem(:,5);

% plot results
hold on;
qmeshcut(cfg.elem(cfg.elemprop>0,1:4),cfg.node,log10(layercw),'y=30','linestyle','none'); view(3)
qmeshcut(cfg.elem(cfg.elemprop>0,1:4),cfg.node,log10(layercw),'z=27','linestyle','none');
box on;
axis equal
view(-56, 22);
