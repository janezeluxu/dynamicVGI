%function [error ] = main(xPoints, yPoints,x1,y1,x2,y2,OrderList,force)
function [error] = main()
%This program solve the 2D diffution equation for variable or unifrom order
%simplex elements, the mesh is structured triangles.

%close all
%clear 
%clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nsd; 
global stabflag;
global itau; 
global kappa; 
global force;
global a;
global TotalDOF;
global BC; 
global minX;
global maxX;
global minY;
global maxY;
global tol; 
global h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [1,0;0,1]*1e-1;
force = 0;
a = [0.5 0];

% kappa = [1,0;0,1];
% force = 0;
% a = [0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 1; 
xPoints = 11;   % number of points in x
yPoints = 11;   % number of points in x
xDomain = [0,1];   % limit of domain in x
yDomain = [0,1];   % limit of domain in x

%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [0];
y1 = [0];
x2 = [1];
y2 = [1];
OrderList = [1];    % define p-order(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-15; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCtop = 0; 
BCbot = 0; 
BCleft = 0;
BCright = 1; 
BC = [BCtop BCbot BCleft BCright]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the mesh Peclet number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nsd==1)
    h = range(xDomain)/(xPoints-1);
    alpha = h*a(1)/(2*kappa(1,1));     
else
    h = sqrt(0.5*range(xDomain)/(xPoints-1) * range(yDomain)/(yPoints-1));
    alpha = h*norm(a)/(2*norm(kappa)); 
end

fprintf('Mesh Peclet Number: %f\n',alpha); 

nInt = 10;

pV = [2];
pH = [1];
pT = [2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mesh data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if force ==7
    [TotalDOF,vertexData,IBC,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,pV,pH,pT); 
else
    [TotalDOF,vertexData,IBC,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList); 
end
u = zeros(TotalDOF,1);    % Step n solution

x = zeros(length(u),1);
y = zeros(length(u),1);
for i = 1:length(vertexData)
    x(i) = vertexData{i,2}(1);
    y(i) = vertexData{i,2}(2);
end
xu = unique(x); 
yu = unique(y); 


tic
%call solveU function, u is the solution field
[u] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,iper,OrderList);

% if force ==7
%     %[u] = solveU(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,kappa,force,pV,pH,pT);
%     [u] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,iper,pV,pH,pT);
%      GridSize = (xPoints-1)*(yPoints-1)*2;
%     %[Error] = errorEstimate(u,xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,force,GridSize,nInt,pV,pH,pT);
%     [Error] = errorEstimate(u,vertexData,IENstruct,GridSize,nInt,pV,pH,pT);
%     error = Error^0.5;
% else
%     [u] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,iper,OrderList);
%      GridSize = (xPoints-1)*(yPoints-1)*2;
%      [Error] = errorEstimate(u,vertexData,IENstruct,GridSize,nInt,OrderList);
%      error = Error^0.5;
%     %mesh = geo_Mesh(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList);
% end

toc


 GridSize = (xPoints-1)*(yPoints-1)*2;
 [Error] = errorEstimate(u,vertexData,IENstruct,GridSize,nInt,OrderList);
 error = Error^0.5;

% %plot the solution using 1D line plot (only vertex modes) or 2D surface plot
figure(1)
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
    c = colorbar;
    title(c,'u','FontSize',25)
    xlabel('x');
    ylabel('y'); 
    zlabel('u'); 
    axis equal
    xlim([minX maxX]);
    ylim([minY maxY]);
    set(gca,'FontSize',15)
    ylabel('y','FontSize',20)
    xlabel('x','FontSize',20)
    zlabel('u','FontSize',20)
    title(horzcat('Pe_{cell} = ',num2str(alpha)),'FontSize',30); 

 [x,y,value] = DiagnalPlot(u,xPoints,vertexData);

  for i = 1:length(x)
     s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
  end
  figure(3)
  plot(s,value,'ro')
  hold on
  uex = Exact(x,y,force,a,kappa);
  value;
  DiagnalError = uex'-value
  plot(s,uex,'b-')
  
  figure(4)
  plot(s,DiagnalError,'b-')
  
end

