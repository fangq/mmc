
[node,face,elem]=meshabox([0 0 0],[60 60 20],2, 2);

data1=load('replay1.dat');
data2=load('replay2.dat');
data3=load('replay3.dat');
data4=load('replay4.dat');
data4_1=data4(1:end/2,:);
data4_2=data4(end/2+1:end,:);

figure(1);
plotmesh([node,log(data1(1:end-7,2))],elem);

figure(2);
plotmesh([node,log(data2(1:end-7,2))],elem);

figure(3);
plotmesh([node,log(data3(1:end-7,2))],elem);

figure(4);
plotmesh([node,log(data4_1(1:end-7,2))],elem);

figure(5);
plotmesh([node,log(data4_2(1:end-7,2))],elem);
