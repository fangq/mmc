function face=readnmrface(filename)
fid=fopen(filename,'rt');
[header,c]=fscanf(fid,'%d',4);
face=fscanf(fid,'%d',[6 header(4)]);
fclose(fid);

face=face';
face=face(:,2:4);
