function elem=readmmcelem(filename)
fid=fopen(filename,'rt');
[header,c]=fscanf(fid,'%d',2);
elem=fscanf(fid,'%d',[6 header(2)]);
fclose(fid);

elem=elem';
elem=elem(:,2:5);
