function [u,c] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,iper,varargin)
%solve Kglobal*u = fGlobal matrix using gmres

global TotalDOF; 
global tol; 

if length(varargin) ==1
    OrderList = varargin{1};
    P = OrderList;
elseif length(varargin) ==3
    PV = varargin{1};
    PH = varargin{2};
    PT = varargin{3};
    P = [PV,PH,PT];
    P = unique(P);
end

[kGlobal, fGlobal] = globalKF(P,IENstruct,IBC,BCval,...
                              vertexData,xPoints,yPoints,iper);

c = condest(kGlobal)  ;                   
u = gmres(kGlobal,fGlobal,50,tol,TotalDOF) ;  
%u = kGlobal\fGlobal;
end