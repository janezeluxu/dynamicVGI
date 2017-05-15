function [] = testLeastSquare(x,y)
x = [1,2,3];
y = [1.2,1,1.1];
[a,b] = LeastSquare(x,y);
[a1,a2] = MatlabLS(x,y);
ae = a-a1;
be = b-a2;

if  (abs(ae) < 1E-14)&&(abs(be) < 1E-14)
   disp( 'test leatSquare Successed!!!!')
else
   disp( 'test leatSquare Failed!!!!')
end
end

function [a1,a2] = MatlabLS(x,u)
node = 2;
uh = [u(node-1);u(node);u(node+1)];
CordMatrix = [1,x(node-1);1,x(node);1,x(node+1)];

A = CordMatrix\uh;
a1 = A(1);
a2 = A(2);
end

function [a,b] = LeastSquare(x,y)
n = length(x);
xbar = sum(x)/n;
ybar = sum(y)/n;
xs = sum(x.^2);
xiyi = sum(x.*y);
a = (ybar*xs-xbar*xiyi)/(xs-n*(xbar^2));
b = (xiyi-n*xbar*ybar)/(xs-n*xbar^2);
end