function detw=mmcdetweight(detpt,prop,w0)
medianum=size(prop,1);
if(medianum<=1)
    error('empty property list');
end
detpt=detpt';
detw=w0';
for i=2:medianum
    detw=detw.*exp(-prop(i,1)*detpt(:,i+medianum-1));
end
end
