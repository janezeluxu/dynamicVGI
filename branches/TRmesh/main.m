%function [error ] = main(xPoints, yPoints,x1,y1,x2,y2,OrderList)
function [] = main(maxIter)
%This program solve the 2D diffution equation for variable or unifrom order
%simplex elements, the mesh is structured triangles.
maxIter = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nsd; 
global stabflag;
global itau; 
global kappa; 
global direction
global force
global BCType;
global a;
global TotalDOF;
global BC; 
global minX;
global maxX;
global minY;
global maxY;
global tol; 
global h;
global Grid_size;
global node_Size;
global gamma;
global IntByPart;
global StepFunc;
global c2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [1,0;0,1]*1e-5;
force = 0;
a = [1 1];

n1 = [0;-1];
n2 = [0;1];
n3 = [-1;0];
n4 = [1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 2; 
xPoints = 41;   % number of points in x
yPoints = 41;   % number of points in x
xDomain = [0,1];   % limit of domain in x
yDomain = [0,1];   % limit of domain in x

%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [0];
y1 = [0];
x2 = [1];
y2 = [1];
OrderList = [1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-9; 
uTol = 1e-3;
%c1Tol = 1e-10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCbot = 0;
BCtop = 0;  
BCleft = 'StepFunc';
BCright = 0;
StepFunc = [0,1,0.3];
BC = {BCbot BCtop BCleft BCright}; 
direction = 'X';
%BCType ={'PeriodicY','','Strong','Strong'};
%BCType = {'Outlet','Inlet','Outlet','Inlet'};
%BCType = {'Outlet','Strong','Strong','Inlet'};
%BCType = {'Strong','Strong','Strong','Strong'};
%BCType = {'flux','Strong','flux','Strong'};
BCType = {'Strong','flux','Strong','flux'};
%BCType = {'flux','flux','flux','Strong'};
IntByPart = 'no';
gamma = 1;
global CbI;
CbI = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 3;  
Grid_size = (xPoints-1)*(yPoints-1)*2;
node_Size = xPoints*yPoints;
nInt = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the mesh Peclet number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nsd==1)
    h = range(xDomain)/(xPoints-1);
    alpha = h*norm(a(1))/(2*kappa(1,1));     
else
    h = sqrt(0.5*range(xDomain)/(xPoints-1) * range(yDomain)/(yPoints-1));
    alpha = h*norm(a)/(2*norm(kappa));
end
%h
fprintf('Mesh Peclet Number: %f\n',alpha); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mesh data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[TotalDOF,vertexData,edgeData,faceData,IBC,uHBCE,uHBCV,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList); 
 
u = zeros(TotalDOF,1);
x = zeros(length(u),1);
y = zeros(length(u),1);
for i = 1:xPoints*xPoints
    x(i) = vertexData{i,2}(1);
    y(i) = vertexData{i,2}(2);
end

xu = unique(x); 
yu = unique(y); 


%call solveU function, u is the solution field
c1 = zeros(node_Size,1);
dt = 3/7;
c2 = 1;
[u] = solveU(vertexData,IBC,BCval,IENstruct,iper,c1,edgeData,n1,n2,n3,n4,OrderList);
c1Start = solvec1(OrderList,xPoints,yPoints,vertexData,faceData,IBC,uHBCE,uHBCV,1,IENstruct,u);
upL2normStart = upnorm(u,c1Start,vertexData,IENstruct,1,OrderList);

L1 = zeros(node_Size,1);
L2 = zeros(node_Size,1);
M1 = zeros(node_Size,1);
M2 = zeros(node_Size,1);
if itau == 3
    for iteration = 1:maxIter
        iteration
        [u] = solveU(vertexData,IBC,BCval,IENstruct,iper,c1,edgeData,n1,n2,n3,n4,OrderList);
        [c1new,L1new,L2new,M1new,M2new] = solvec1(OrderList,xPoints,yPoints,vertexData,faceData,IBC,uHBCE,uHBCV,1,IENstruct,u);
        
        upL2norm1 = upnorm(u,c1,vertexData,IENstruct,1,OrderList);
        [c1,L1,L2,M1,M2] = update_c1(c1,L1,L2, M1,M2, L1new,L2new,M1new,M2new, dt,xPoints,yPoints,vertexData,IBC,'pathLinec1');
        %c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
        
        [u2] = solveU(vertexData,IBC,BCval,IENstruct,iper,c1,edgeData,n1,n2,n3,n4,OrderList);
        upL2norm2 = upnorm(u2,c1,vertexData,IENstruct,1,OrderList);
        normu = norm(u-u2)/norm(u)
        erru = norm(upL2norm2- upL2norm1)/norm(upL2normStart)
        
        if  ( normu<uTol) %erru <uTol ||
            break
        end
    end
end
 [Error] = errorEstimate(u,vertexData,IENstruct,nInt,OrderList);
 error = Error^0.5;

%plot the solution using 1D line plot (only vertex modes) or 2D surface plot
figure(2)
count = 1; 
xu = unique(x); 
yu = unique(y);  
zz = zeros(length(xu),length(yu)); 
for j = 1:length(yu)
    for i = 1:length(xu)
        zz(j,i) = u(count); 
        count = count+1; 
    end
end

grid on

minX = min(xDomain); 
maxX = max(xDomain); 
minY = min(yDomain); 
maxY = max(yDomain); 

surf(xu,yu,zz)
view([-40 30]);
axis([minX maxX minY maxY -0.2 1.2])
[x,y,value] = DiagnalPlot(u,xPoints,vertexData);

 for i = 1:length(x)
     s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
 end
 figure(10)
 plot(s,value,'r*')
 
 hold on
 uex = Exact(x,y,force,a,kappa);
 value;
 DiagnalError = uex'-value;
 %plot(s,uex,'b-')
  
  
end

