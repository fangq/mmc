
[node,face,elem]=meshabox([0 0 0],[60 60 20],2, 2);

data1=load('replay1.dat');
data2=load('replay2.dat');
data3=load('replay3.dat');

figure(1);
plotmesh([node,log(data1(1:end-7,2))],elem);

figure(2);
plotmesh([node,log(data2(1:end-7,2))],elem);

figure(3);
plotmesh([node,log(data3(1:end-7,2))],elem);