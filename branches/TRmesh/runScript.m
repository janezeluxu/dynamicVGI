
global nsd
global direction
global force
global BCType;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 1; 
force = 6;
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1,-0.6];
y1 = [-1,-1];
x2 = [1,0.6];
y2 = [1,1];
OrderList = [3,3];

[error6 ] = main(6,6,x1,y1,x2,y2,OrderList);
fprintf('error: %d\n',error6);
[error11 ] = main(11,11,x1,y1,x2,y2,OrderList);
fprintf('error: %d\n',error11);
% [error21 ] = main(21,21,x1,y1,x2,y2,OrderList);
% fprintf('error: %d\n',error21);
% [error41 ] = main(41,41,x1,y1,x2,y2,OrderList);
% fprintf('error: %d\n',error41);
%[error81 ] = main(81,81,x1,y1,x2,y2,OrderList);
%fprintf('error: %d\n',error81);
