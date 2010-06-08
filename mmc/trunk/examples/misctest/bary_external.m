% test bary centric coordinates outside a tetrahedron

node=[0 0 0;1 0 0;0 1 0;0 0 1];
elem=[1 2 3 4];
tetramesh(elem,node);
view([60,30]);
hold on;
b1=[1/3,1/3,1/3,0];
b2=[0 1/3 1/3 1/3];
t=-2:0.1:2;

for i=1:length(t)
 pp=node'*(t(i)*b1+(1-t(i))*b2)';
 plot3(pp(1),pp(2),pp(3),'o');
 axis equal;
 waitforbuttonpress;
end

