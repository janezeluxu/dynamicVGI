function [u] = solveU(vertexData,IBC,BCval,IENstruct,iper,c1,EdgeData,n1,n2,n3,n4,varargin)
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
                              vertexData,iper,c1,EdgeData,n1,n2,n3,n4);
c = condest(kGlobal);                   
%u = gmres(kGlobal,fGlobal,50,tol,TotalDOF) ;  
u = kGlobal\fGlobal;
end