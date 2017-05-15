function [error] = ErrorEstimate(u,nInt,Grid_Size,a,kappa)
global h

nElement = Grid_Size-1;
error = 0;
[xi,weight] = GaussQuad(nInt, 1);
%nElement = Grid_Size-1;
for ele = 1:nElement
    uE = [u(ele),u(ele+1)];
    len = length(xi);
    ShapeFunc = zeros(2,len);
    %uA = zeros(1,len);
    for n = 1:len
        ShapeFunc(:,n) = [1-xi(n);xi(n)];
    end
    uH = uE*ShapeFunc;
    
    x =  (ele-1+xi).*h(ele);
    uA = Exact(x,a,kappa);
    error = ((uH-uA').^2*weight*h(ele))+error;
    %error = error.^0.5
    
end
error = error^0.5;
end

