function [] = testError()
global h 
Grid_Size = 9;
for i = 1:Grid_Size-1
    h(i) = 1/(Grid_Size-1);
end
u = zeros(Grid_Size,1);
[Error] = ErrorEstimate(u,60,Grid_Size,-1,1e-3);
ErrorA = analyticError(Grid_Size,-1,1e-3);
double(ErrorA)

errot  = Error - double(ErrorA);
L2 = (sum(sum(errot.^2)))^0.5;
if  (L2 < 1E-15)
   disp( 'test error estimate Successed!!!!')
else
   disp( 'test error estimate Failed!!!!')
end

end
function IntEx = analyticError(Grid_Size,a,kappa)
IntEx = 0;

nElement = Grid_Size-1;
h = zeros(1,nElement);
for i = 1:nElement
    h(i) = (1)/(nElement);
end

for ele = 1:nElement
    syms x
    x1 = (ele-1)*h(ele);
    x2 = (ele)*h(ele);
    uex = Exact(x,a,kappa);
    IntEx = int(uex^2,x,x1,x2)+IntEx;
end
IntEx = IntEx^0.5;
end