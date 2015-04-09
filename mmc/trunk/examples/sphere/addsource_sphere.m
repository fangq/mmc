[node,face,elem]=meshasphere([0 0 0],24,1,2);
elem(:,5)=1;
 
source=double([-5 -5;-5 5;5 5;5 -5]);
source(:,3)=25;
[newnode,newelem]=mmcaddsrc(node,elem,source);
plotmesh(newnode,newelem,'x>-5');

savemmcmesh('sphere',newnode,newelem);
