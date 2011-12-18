addpath('../../matlab/')
chpt=[1e4 1e5 1e6 1e7];
no=readmmcnode('node_cube20.dat');

data=zeros(length(chpt),10,size(no,1));
for i=1:length(chpt)
      for j=1:10 
          dd=load(sprintf('v_%d_%d.dat',j,chpt(i)));
          data(i,j,:)=dd(:,end);
      end
end

data=log10(data);
%chpt=sqrt(chpt);

idx=find(ismember(no,[10 10 0],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'bo-');
hold on;
idx=find(ismember(no,[10 10 5],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'ro-');
idx=find(ismember(no,[10 10 10],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'ko-');
idx=find(ismember(no,[10 10 15],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'co-');
idx=find(ismember(no,[10 5 5],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'bo:');
idx=find(ismember(no,[10 5 10],'rows'));
loglog(chpt,[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))],'ro:');

