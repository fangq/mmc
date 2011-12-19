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

%data=log10(data);
chpt=sqrt(chpt);

figure; hold on;
sty={'bo-','ro-','ko-','co-','bo:','ro:'};
pos=[10 10 0;10 10 5;10 10 10;10 10 15;10 5 5;10 5 10];

% plot trial numbers with Coefficient of Variation (100*std/mean)

for i=1:size(pos,1)
	idx=find(ismember(no,pos(i,:),'rows'));
	plot(chpt,100*[std(data(1,:,idx)),std(data(2,:,idx)),std(data(3,:,idx)),std(data(4,:,idx))]./chpt ... % standard error
                    ./[mean(data(1,:,idx)),mean(data(2,:,idx)),mean(data(3,:,idx)),mean(data(4,:,idx))],sty{i});
end
hold off;
