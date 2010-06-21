function node=readmmcnode(filename)

fid=fopen(filename,'rt');
hd=fscanf(fid,'%d',2);
node=fscanf(fid,'%d %e %e %e\n',hd(2)*4);
fclose(fid);

node=reshape(node,[4,hd(2)])';
node=node(:,2:4);
