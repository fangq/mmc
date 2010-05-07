dim=61;
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
dist=(xi-31).^2+(yi-31).^2+(zi-31).^2;

v=zeros(size(xi));
v(dist<400)=1;
v=v+1;

fid=fopen('spherebox.bin','wb');
aa=fwrite(fid,v,'uchar');
fclose(fid);
