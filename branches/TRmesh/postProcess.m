function [L2Error] = postProcess(u,xPoints,yPoints,xDomain,yDomain,...
                x1,y1,x2,y2,OrderList,xLine,yLine,nPoints)

[x,y,value] = Reconstruct(u,xPoints,yPoints,xDomain,yDomain,...
                x1,y1,x2,y2,OrderList,xLine,yLine,nPoints);
 for i = 1:length(x)
    s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
 end
 figure(2)
 plot(s,value,'bo')
 hold on
 [uex] = SteadyHeatAnalytical(x,y);
 
 %error estimate 
 error = uex-value;
 L2Error = (sum(error.^2))^0.5;
 
 %plot analytical solution
 x = sort(x);
 y = sort(y,'descend');
 for i = 1:length(x)
    s(i) = ((x(i) - x(1))^2 + (y(i) - y(1))^2)^0.5;
 end
 [uex] = SteadyHeatAnalytical(x,y);
  plot(s,uex,'b-')

  figure(3)
  hold on
  plot(s,error,'go-')
  ylim([-0.0001 0.0002])
  hold on
  grid on
  
end