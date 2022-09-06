% read output
fid=fopen('optix.bin');
output=fread(fid,'float64');

% retrieve time-resolved results(10 time gates)
res=reshape(output,[61,61,61,10]);

% convert to cw solution and visualize
cw=sum(res,4);
mcxplotvol(log(cw));