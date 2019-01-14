
% [node,face,elem]=meshabox([0 0 0],[60 60 20],2, 2);
nTG = 50;
npat = 5;

rawdata=load('sharing.dat');
data=reshape(rawdata(:,end),npat,[],nTG);
cwdata=squeeze(sum(data,3));
cwdata=cwdata(:,1:end-3);

%% plot results

figure();

subplot(1,3,1); title('pattern 1');
plotmesh([node,(cwdata(1,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,2); title('pattern 2');
plotmesh([node,(cwdata(2,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,3); title('pattern 3');
plotmesh([node,(cwdata(3,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,5); title('pattern 4');
plotmesh([node,(cwdata(4,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])

subplot(2,3,6); title('pattern 5');
plotmesh([node,(cwdata(5,:))'],elem,'linestyle','none');
shading interp
axis off
view([0,0,-1])