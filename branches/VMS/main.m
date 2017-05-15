% function [ ] = main(  )
%This program solve the 2D diffution equation for variable or unifrom order
%simplex elements, the mesh is structured triangles.

close all
clear 
clc

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
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 1; 
xPoints = 4;   % number of points in x
yPoints = 4;   % number of points in x
xDomain = [0,1];   % limit of domain in x
yDomain = [0,1];   % limit of domain in x
%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [0];
y1 = [0];
x2 = [1];
y2 = [1];
OrderList = [1];    % define p-order(s)
minX = min(xDomain); 
maxX = max(xDomain); 
minY = min(yDomain); 
maxY = max(yDomain); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [.1,0;0,.1];
force = 0;
a = [.5 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCtop = 0; 
BCbot = 0; 
BCleft = 0;
BCright = 1; 
BC = [BCtop BCbot BCleft BCright]; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-15; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 1;       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mesh data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[TotalDOF,vertexData,IBC,BCval,IENstruct,iper] = buildMeshStruct(xPoints,yPoints,xDomain,yDomain,x1,y1,x2,y2,OrderList); 

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

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call solveU function, u is the solution field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(TotalDOF,1);    % Step n solution

x = zeros(length(u),1);
y = zeros(length(u),1);
for i = 1:length(vertexData)
    x(i) = vertexData{i,2}(1);
    y(i) = vertexData{i,2}(2);
end
xu = unique(x); 
yu = unique(y); 
    
[u] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,OrderList,iper);
    
% Plot
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

if(nsd==1)
    % 1D Plots   
    x1 = zeros(xPoints,1);
    u1 = zeros(xPoints,1); 
    u2 = zeros(xPoints,1); 
    count1 = 1;     
    count2 = 1;     
    for i = 1:length(y)
        if(y(i)==min(yDomain))
            x1(count1) = x(i);
            u1(count1) = u(i); 
            count1 = count1+1; 
        end   
    end
    
    % Exact solutiion
    L = range(xDomain); 
    xt = (x1-min(x1))/L;     
    Pe = a(1)*L/kappa(1,1); 
    ue = BCleft*((1-exp(Pe*(xt-1)))/(1-exp(-Pe))) + BCright*((exp(Pe*(xt-1))-exp(-Pe))/(1-exp(-Pe))) + (force*L/a(1))*(xt + (exp(-Pe)-exp(Pe*(xt-1)))/(1-exp(-Pe)));

    error = norm(u1-ue)/norm(ue); 
    fprintf('Error: %f\n',error); 
    
    plot(x1,ue,'kx','LineWidth',2.5,'MarkerSize',10)
    hold on
    plot(x1,u1,'b','LineWidth',2.5,'MarkerSize',10)
    legend('Exact','VMS','Location','NorthWest')
    xlabel('x','FontSize',15);
    ylabel('u','FontSize',15); 
    xlim([minX maxX]);
%     ylim([0 1]);
    title(horzcat('Pe = ',num2str(Pe),', Pe_{cell} = ',num2str(alpha),', N_{el} = ',num2str(xPoints-1)),'FontSize',25); 
    ax = gca; 
    set(ax,'XTick',min(xDomain):range(xDomain)/(xPoints-1):max(xDomain));
    grid on
    set(gca,'FontSize',15)
    
    figure(2)
    plot(x1,abs(ue-u1),'b-x','LineWidth',2.5,'MarkerSize',10)
    grid on
%     ylim([1e-16 1])
    xlim([0 1])
    xlabel('x','FontSize',15)
    ylabel('|u-u_{exact}|','FontSize',15)
    title(horzcat('Pe_{cell} = ',num2str(alpha),', N_{el,x} = ',num2str(xPoints-1),', N_{el,y} = ',num2str(yPoints-1)),'FontSize',25); 
    set(gca,'FontSize',15)
else
    %2D Plots    
    figure(2)
    subplot(1,2,1)
    contourf(xu,yu,zz,'LineStyle','None')
    c = colorbar;
    title(c,'u','FontSize',25)  
    caxis([-.1 1.3])
    xlabel('x');
    ylabel('y'); 
    zlabel('u'); 
    axis equal
    xlim([minX maxX]);
    ylim([minY maxY]);
    title(horzcat('Pe = ',num2str(alpha))); 
    ax = gca; 
    set(ax,'XTick',min(xDomain):range(xDomain)/(xPoints-1):max(xDomain));
    set(ax,'YTick',min(yDomain):range(yDomain)/(yPoints-1):max(yDomain));
    grid on

    subplot(1,2,2)
    surf(xu,yu,zz)
    c = colorbar;
    title(c,'u','FontSize',25)
    caxis([-.1 1.3])
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
end   
    
toc

% end

