function [u] = SteadyHeatAnalytical(x,y)
% xs = x(1);
% xt = x(2);
% ys = y(1);
% yt = y(2);
% k = (ys-yt)/(xs-xt);
% b = ys-k*xs;
% dx = (2)/nPoints;
% x = -1:dx:1;
% y = k*x+b;
PI = 3.1415926;
 sum = 0;
 for k = 1:2:11
     devide = k^3*sinh(k*PI);
     times = sinh(k*PI*(1+y)./2)+sinh(k*PI*(1-y)./2);
     xPart = sin(k*PI*(1+x)./2);
     a = (xPart./devide).*(times);
     sum = sum+a;
 end
 %sum
  u  = (1-x.^2)./2 - (16/(PI^3))*(sum);
  
end