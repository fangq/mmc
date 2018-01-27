function detidx = mmcdetidx(detpos, pos0, detparam1, detparam2)

% During wide-field detection, calculate the index of detected photons
% on a 2D detection array.
% detp stores all the information of detected photons
% detpos specifies the starting point of the 2D detection plane
% detparam1&2 specifies the size and resolution of the detection plane

x0 = pos0(1,1);
y0 = pos0(1,2);

xrange = detparam1(1) + detparam2(1);
yrange = detparam1(1) + detparam2(1);

xsize = detparam1(4);
ysize = detparam2(4);

xindex = floor((detpos(:,1)-x0) / xrange * xsize);
yindex = floor((detpos(:,2)-y0) / yrange * ysize);

assert(isempty(find(xindex<0, 1) & find(xindex>=xsize, 1)), "photon location not within detection plane");
assert(isempty(find(yindex<0, 1) & find(yindex>=ysize, 1)), "photon location not within detection plane");

detidx = single(yindex * ysize + xindex);

end

