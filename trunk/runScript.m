%close all
clear 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [1,0;0,1]*1e-5;
force = 0;
a = [-0.5 -1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 2; 
xPoints = 11;   % number of points in x
yPoints = 11;   % number of points in x
xDomain = [-1,1];   % limit of domain in x
yDomain = [-1,1];   % limit of domain in x

%define the order list of domain list from [x1(i),y1(i)] to [x2(i),y2(i)]
x1 = [-1];
y1 = [-1];
x2 = [1];
y2 = [1];
OrderList = [1];    % define p-order(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-12; 

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
itau = 2; 

[error] = main(xDomain, yDomain,xPoints, yPoints,x1,y1,x2,y2,OrderList);
