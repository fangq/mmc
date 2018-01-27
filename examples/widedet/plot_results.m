
nx = 40; ny=40;
nTG = 50;

fd=fopen('widedet.img','rb','ieee-le');
data = fread(fd,nx*nx*nTG,'float');
data=reshape(data,[nx ny nTG]);
fclose(fd);

datacw=sum(data,3);
figure(1);
imagesc(datacw');
axis equal;
xlim([0.5 nx+0.5]);
xlabel('offsetx (mm)');
ylabel('offsety (mm)');
