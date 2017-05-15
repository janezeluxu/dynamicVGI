function [u] = solveU(xPoints,yPoints,vertexData,IBC,BCval,IENstruct,OrderList,iper)

global TotalDOF; 
global tol; 

[kGlobal, fGlobal] = globalKF(OrderList,IENstruct,IBC,BCval,...
                              vertexData,xPoints,yPoints,iper)
 
u = gmres(kGlobal,fGlobal,50,tol,TotalDOF) ;  

end